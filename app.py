"""
Primerly - PCR Primer Design Tool
A Biovagon Initiative
Copyright (c) 2025 Jeevitha C M, Founder of Biovagon

A beginner-friendly web application for designing PCR primers using primer3-py and Biopython.
Provides an intuitive interface for both standard PCR and qPCR primer design with comprehensive quality analysis.

Version: 1.0.0
License: MIT
Author: Jeevitha C M (Founder of Biovagon)
Website: https://www.biovagon.org/
"""

from flask import Flask, render_template, request, jsonify, send_file, redirect, session
from Bio import Entrez
from Bio.Seq import Seq
import primer3
import csv
import io
import re
import uuid
from datetime import datetime
from collections import defaultdict

# Application metadata
__version__ = "1.0.0"
__author__ = "Jeevitha C M"
__organization__ = "Biovagon"
__license__ = "MIT"
__copyright__ = "Copyright (c) 2025 Jeevitha C M, Founder of Biovagon"
__website__ = "https://www.biovagon.org/"

app = Flask(__name__)
app.secret_key = 'your-secret-key-here'

# In-memory storage for results (use database for production)
result_store = {}

# Set email for NCBI Entrez (required)
Entrez.email = "jeevithacm21@gmail.com"

def fetch_sequence_from_ncbi(accession):
    """Fetch DNA sequence from NCBI using accession number"""
    try:
        # Search for the accession number
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        
        if not record or record.startswith("Error"):
            return None, "Invalid accession number or sequence not found"
        
        # Parse FASTA format
        lines = record.strip().split('\n')
        if len(lines) < 2:
            return None, "Invalid sequence format"
        
        # Extract sequence (join all lines except the first which is the header)
        sequence = ''.join(lines[1:]).upper()
        
        # Validate sequence contains only valid DNA characters
        if not re.match(r'^[ATCGN]+$', sequence):
            return None, "Sequence contains invalid characters"
        
        return sequence, None
    except Exception as e:
        return None, f"Error fetching sequence: {str(e)}"

def calculate_primer_score(forward_primer, reverse_primer, product_size, experiment_type='standard_pcr'):
    """Calculate a score (0-100) for primer pair quality with specific explanations"""
    score = 0
    explanations = []
    
    # Tm between 58-62°C (20 pts) - tighter for qPCR
    forward_tm = primer3.calc_tm(forward_primer)
    reverse_tm = primer3.calc_tm(reverse_primer)
    
    if experiment_type == 'qpcr':
        # qPCR: tighter Tm range
        if 58 <= forward_tm <= 62 and 58 <= reverse_tm <= 62:
            score += 20
            explanations.append(f"Perfect qPCR Tm: F={forward_tm}°C, R={reverse_tm}°C (optimal 58-62°C)")
        elif 57 <= forward_tm <= 63 and 57 <= reverse_tm <= 63:
            score += 10
            explanations.append(f"Good qPCR Tm: F={forward_tm}°C, R={reverse_tm}°C (acceptable 57-63°C)")
        else:
            explanations.append(f"Tm outside qPCR range: F={forward_tm}°C, R={reverse_tm}°C (need 57-63°C)")
    else:
        # Standard PCR: more flexible
        if 58 <= forward_tm <= 62 and 58 <= reverse_tm <= 62:
            score += 20
            explanations.append(f"Optimal Tm: F={forward_tm}°C, R={reverse_tm}°C (perfect 58-62°C)")
        elif 55 <= forward_tm <= 65 and 55 <= reverse_tm <= 65:
            score += 10
            explanations.append(f"Acceptable Tm: F={forward_tm}°C, R={reverse_tm}°C (good 55-65°C)")
        else:
            explanations.append(f"Tm outside range: F={forward_tm}°C, R={reverse_tm}°C (need 55-65°C)")
    
    # GC content between 40-60% (15 pts)
    forward_gc = (forward_primer.count('G') + forward_primer.count('C')) / len(forward_primer) * 100
    reverse_gc = (reverse_primer.count('G') + reverse_primer.count('C')) / len(reverse_primer) * 100
    
    if 40 <= forward_gc <= 60 and 40 <= reverse_gc <= 60:
        score += 15
        explanations.append(f"Balanced GC: F={forward_gc:.1f}%, R={reverse_gc:.1f}% (optimal 40-60%)")
    elif 35 <= forward_gc <= 65 and 35 <= reverse_gc <= 65:
        score += 7
        explanations.append(f"Moderate GC: F={forward_gc:.1f}%, R={reverse_gc:.1f}% (acceptable 35-65%)")
    else:
        explanations.append(f"GC content issue: F={forward_gc:.1f}%, R={reverse_gc:.1f}% (need 35-65%)")
    
    # Primer length between 18-25 bp (10 pts)
    forward_len = len(forward_primer)
    reverse_len = len(reverse_primer)
    if 18 <= forward_len <= 25 and 18 <= reverse_len <= 25:
        score += 10
        explanations.append(f"Optimal length: F={forward_len}bp, R={reverse_len}bp (perfect 18-25bp)")
    elif 16 <= forward_len <= 30 and 16 <= reverse_len <= 30:
        score += 5
        explanations.append(f"Good length: F={forward_len}bp, R={reverse_len}bp (acceptable 16-30bp)")
    else:
        explanations.append(f"Length issue: F={forward_len}bp, R={reverse_len}bp (need 16-30bp)")
    
    # GC clamp at 3' end (10 pts)
    forward_gc_clamp = forward_primer[-1] in 'GC'
    reverse_gc_clamp = reverse_primer[-1] in 'GC'
    if forward_gc_clamp and reverse_gc_clamp:
        score += 10
        explanations.append(f"Strong GC clamps: F ends with '{forward_primer[-1]}', R ends with '{reverse_primer[-1]}'")
    elif forward_gc_clamp or reverse_gc_clamp:
        score += 5
        explanations.append(f"Partial GC clamp: F ends with '{forward_primer[-1]}', R ends with '{reverse_primer[-1]}'")
    else:
        explanations.append(f"No GC clamp: F ends with '{forward_primer[-1]}', R ends with '{reverse_primer[-1]}'")
    
    # Hairpin structure ∆G > -2 kcal/mol (15 pts)
    forward_hairpin = primer3.calc_hairpin(forward_primer)
    reverse_hairpin = primer3.calc_hairpin(reverse_primer)
    
    if forward_hairpin.dg > -2 and reverse_hairpin.dg > -2:
        score += 15
        explanations.append(f"Clean structure: F hairpin ∆G={forward_hairpin.dg:.1f}, R hairpin ∆G={reverse_hairpin.dg:.1f}")
    elif forward_hairpin.dg > -3 and reverse_hairpin.dg > -3:
        score += 7
        explanations.append(f"Minor hairpins: F hairpin ∆G={forward_hairpin.dg:.1f}, R hairpin ∆G={reverse_hairpin.dg:.1f}")
    else:
        explanations.append(f"Hairpin risk: F hairpin ∆G={forward_hairpin.dg:.1f}, R hairpin ∆G={reverse_hairpin.dg:.1f}")
    
    # No 3' self-dimers (10 pts)
    forward_dimer = primer3.calc_homodimer(forward_primer)
    reverse_dimer = primer3.calc_homodimer(reverse_primer)
    
    if forward_dimer.dg > -2 and reverse_dimer.dg > -2:
        score += 10
        explanations.append(f"No self-dimers: F dimer ∆G={forward_dimer.dg:.1f}, R dimer ∆G={reverse_dimer.dg:.1f}")
    elif forward_dimer.dg > -3 and reverse_dimer.dg > -3:
        score += 5
        explanations.append(f"Minor self-dimers: F dimer ∆G={forward_dimer.dg:.1f}, R dimer ∆G={reverse_dimer.dg:.1f}")
    else:
        explanations.append(f"Self-dimer risk: F dimer ∆G={forward_dimer.dg:.1f}, R dimer ∆G={reverse_dimer.dg:.1f}")
    
    # Product size scoring based on experiment type (10 pts)
    if experiment_type == 'qpcr':
        if 80 <= product_size <= 200:
            score += 10
            explanations.append(f"Perfect qPCR size: {product_size}bp (optimal 80-200bp)")
        elif 60 <= product_size <= 250:
            score += 5
            explanations.append(f"Good qPCR size: {product_size}bp (acceptable 60-250bp)")
        else:
            explanations.append(f"qPCR size issue: {product_size}bp (need 60-250bp)")
    else:
        if 200 <= product_size <= 1000:
            score += 10
            explanations.append(f"Optimal PCR size: {product_size}bp (perfect 200-1000bp)")
        elif 100 <= product_size <= 1200:
            score += 5
            explanations.append(f"Good PCR size: {product_size}bp (acceptable 100-1200bp)")
        else:
            explanations.append(f"PCR size issue: {product_size}bp (need 100-1200bp)")
    
    # Bonus for clean secondary structure (10 pts)
    if score >= 70:
        score += 10
        explanations.append("Excellent overall quality - highly recommended")
    
    return min(score, 100), explanations

def get_score_label(score):
    """Get color label for score"""
    if score >= 80:
        return "Excellent"
    elif score >= 60:
        return "Good"
    else:
        return "Poor"

def visualize_dimer(seq1, seq2, dimer_type="hetero"):
    """
    Visualize dimer formation between two DNA sequences using a sliding window alignment.
    Args:
        seq1: First DNA sequence
        seq2: Second DNA sequence (or same as seq1 for self-dimer)
        dimer_type: "self" or "hetero"
    Returns:
        Dictionary with visualization data
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    best_score = 0
    best_offset = 0
    best_direction = 'left'  # 'left' means seq2 is left-padded, 'right' means seq1 is left-padded
    best_bonds = ''
    best_seq1_aligned = ''
    best_seq2_aligned = ''
    n1 = len(seq1)
    n2 = len(seq2)
    # Slide seq2 left over seq1
    for offset in range(-n2+1, n1):
        s1 = []
        s2 = []
        bonds = []
        score = 0
        for i in range(max(n1, n2 + abs(offset))):
            c1 = seq1[i] if 0 <= i < n1 else ' '
            c2 = seq2[i - offset] if 0 <= i - offset < n2 else ' '
            s1.append(c1)
            s2.append(c2)
            if c1 in complement and complement[c1] == c2:
                bonds.append('|')
                score += 1
            else:
                bonds.append(' ')
        if score > best_score:
            best_score = score
            best_offset = offset
            best_direction = 'left'
            best_bonds = ''.join(bonds)
            best_seq1_aligned = ''.join(s1)
            best_seq2_aligned = ''.join(s2)
    # If no base pairs found
    if best_score == 0:
        return {
            'has_dimer': False,
            'visualization': "No significant dimer formation detected",
            'delta_g': 0,
            'base_pairs': 0
        }
    delta_g = -2.0 * best_score
    visualization = f"Delta G: {delta_g:.2f} kcal/mol  Base Pairs: {best_score}\n"
    visualization += f"5'  {best_seq1_aligned}\n"
    visualization += f"    {best_bonds}\n"
    visualization += f"3'  {best_seq2_aligned}"
    return {
        'has_dimer': True,
        'visualization': visualization,
        'delta_g': delta_g,
        'base_pairs': best_score
    }

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/fetch-sequence', methods=['POST'])
def fetch_sequence():
    """Fetch sequence from NCBI accession number"""
    data = request.get_json()
    accession = data.get('accession', '').strip()
    
    if not accession:
        return jsonify({'error': 'Please enter an accession number'})
    
    sequence, error = fetch_sequence_from_ncbi(accession)
    
    if error:
        return jsonify({'error': error})
    
    return jsonify({'sequence': sequence})

def get_highlighted_sequence(sequence, primers, primer_index=0, orf_info=None, utr_info=None):
    """Return a list of dicts: {char, color_class} for each base in sequence, highlighting specific primer positions."""
    seq_len = len(sequence)
    highlight = [''] * seq_len
    
    # Define colors for each primer pair
    primer_colors = [
        {'forward': 'bg-primary', 'reverse': 'bg-danger'},      # Blue/Red for primer 1
        {'forward': 'bg-success', 'reverse': 'bg-warning'},     # Green/Orange for primer 2  
        {'forward': 'bg-info', 'reverse': 'bg-secondary'}       # Cyan/Gray for primer 3
    ]
    
    # Highlight ORF region if provided
    if orf_info:
        orf_start = orf_info['start']
        orf_end = orf_info['end']
        for i in range(orf_start, orf_end):
            if 0 <= i < seq_len:
                highlight[i] = 'orf-highlight'  # Use the new CSS class with background color
    
    # Highlight 3' UTR region if provided
    if utr_info:
        utr_start = utr_info['start']
        utr_end = utr_info['end']
        for i in range(utr_start, utr_end):
            if 0 <= i < seq_len:
                highlight[i] = 'utr-highlight'  # Use a new CSS class for UTR highlighting
    
    # For each primer, mark the specific positions calculated by primer3
    for idx, primer in enumerate(primers):
        if idx >= len(primer_colors):
            break  # Only handle up to 3 primers
            
        # Get colors for this primer pair - use primer_index for individual sequences
        if len(primers) == 1:
            # Individual sequence highlighting - use the specified primer_index
            forward_color = primer_colors[primer_index]['forward']
            reverse_color = primer_colors[primer_index]['reverse']
        else:
            # Multiple primers - use the index in the list
            forward_color = primer_colors[idx]['forward']
            reverse_color = primer_colors[idx]['reverse']
        
        # Debug: Print primer information
        print(f"Primer {idx + 1} (index {primer_index}):")
        print(f"  Forward: {primer['forward']} (length: {len(primer['forward'])})")
        print(f"  Reverse: {primer['reverse']} (length: {len(primer['reverse'])})")
        print(f"  Colors: {forward_color} / {reverse_color}")
        
        # Highlight forward primer positions (direct sequence match)
        if 'forward_start' in primer and 'forward_end' in primer:
            forward_start = primer['forward_start']
            forward_end = primer['forward_end']
            print(f"  Forward positions: {forward_start}-{forward_end}")
            
            # Highlight the forward primer positions
            for i in range(forward_start, forward_end):
                if 0 <= i < seq_len:
                    # Combine ORF highlighting with primer highlighting
                    if highlight[i] == 'orf-highlight':
                        highlight[i] = f"{forward_color} orf-highlight"
                    else:
                        highlight[i] = forward_color
                    print(f"    Highlighting position {i}: {sequence[i]} with {highlight[i]}")
        
        # Highlight reverse primer positions (direct sequence match)
        if 'reverse_start' in primer and 'reverse_end' in primer:
            reverse_start = primer['reverse_start']
            reverse_end = primer['reverse_end']
            print(f"  Reverse positions: {reverse_start}-{reverse_end}")
            
            # Highlight the reverse primer positions
            for i in range(reverse_start, reverse_end):
                if 0 <= i < seq_len:
                    # Combine ORF highlighting with primer highlighting
                    if highlight[i] == 'orf-highlight':
                        highlight[i] = f"{reverse_color} orf-highlight"
                    else:
                        highlight[i] = reverse_color
                    print(f"    Highlighting position {i}: {sequence[i]} with {highlight[i]}")
    
    # Convert to list of dicts
    result = []
    for i, char in enumerate(sequence):
        result.append({
            'char': char,
            'color_class': highlight[i] if highlight[i] else ''
        })
    
    return result

@app.route('/result/<result_id>')
def show_result(result_id):
    """Display primer design results using stored data"""
    try:
        # Get results from storage
        data = result_store.get(result_id)
        if not data:
            return "Result not found or expired", 404
        
        print(f"DEBUG: Showing result for ID: {result_id}")
        print(f"DEBUG: Found {len(data['primers'])} primers")
        
        return render_template('results.html', 
                             primers=data['primers'], 
                             sequence_length=data['sequence_length'],
                             experiment_type=data['experiment_type'],
                             sequence=data['sequence'],
                             highlighted_sequence=data['highlighted_sequence'],
                             individual_highlighted_sequences=data['individual_highlighted_sequences'],
                             orf_info=data['orf_info'],
                             utr_info=data['utr_info'],
                             mode=data['mode'])
    except Exception as e:
        print(f"Error in result route: {e}")
        import traceback
        traceback.print_exc()
        return f"Error processing result: {str(e)}", 500

@app.route('/results')
def results():
    """Legacy results route - redirects to new system"""
    return "This route is deprecated. Please use the new result system.", 400

def find_largest_orf(sequence):
    """
    Find the largest Open Reading Frame (ORF) in the given DNA sequence.
    Returns the start position, end position, and the ORF sequence.
    """
    try:
        # Define start and stop codons
        start_codons = ['ATG']  # Only ATG start codon
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        print(f"DEBUG: Looking for ORFs with ATG start codon only")
        
        # Find all possible ORFs in all three reading frames
        orfs = []
        
        for frame in range(3):  # Three reading frames
            frame_sequence = sequence[frame:]
            
            # Find all start codons in this frame
            start_positions = []
            for i in range(0, len(frame_sequence) - 2, 3):
                codon = frame_sequence[i:i+3]
                if codon in start_codons:
                    start_positions.append(i)
            
            print(f"DEBUG: Frame {frame + 1}: Found {len(start_positions)} ATG start codons")
            
            # For each start codon, find the next stop codon
            for start_idx in start_positions:
                start_pos = frame + start_idx  # Convert to original sequence position
                
                # Look for stop codon after this start codon
                for j in range(start_idx + 3, len(frame_sequence) - 2, 3):
                    codon = frame_sequence[j:j+3]
                    if codon in stop_codons:
                        end_pos = frame + j + 3  # Convert to original sequence position
                        orf_length = end_pos - start_pos
                        
                        # Only include ORFs longer than 30 nucleotides (10 codons)
                        if orf_length >= 30:
                            orf_sequence = sequence[start_pos:end_pos]
                            orfs.append({
                                'start': start_pos,
                                'end': end_pos,
                                'length': orf_length,
                                'sequence': orf_sequence,
                                'frame': frame + 1
                            })
                        break
        
        # Also check the reverse complement
        try:
            reverse_complement = str(Seq(sequence).reverse_complement())
            
            for frame in range(3):  # Three reading frames on reverse complement
                frame_sequence = reverse_complement[frame:]
                
                # Find all start codons in this frame
                start_positions = []
                for i in range(0, len(frame_sequence) - 2, 3):
                    codon = frame_sequence[i:i+3]
                    if codon in start_codons:
                        start_positions.append(i)
                
                # For each start codon, find the next stop codon
                for start_idx in start_positions:
                    # Convert reverse complement position to original sequence position
                    start_pos_rc = frame + start_idx
                    start_pos = len(sequence) - start_pos_rc - 3  # Convert to original sequence
                    
                    # Look for stop codon after this start codon
                    for j in range(start_idx + 3, len(frame_sequence) - 2, 3):
                        codon = frame_sequence[j:j+3]
                        if codon in stop_codons:
                            end_pos_rc = frame + j + 3
                            end_pos = len(sequence) - end_pos_rc  # Convert to original sequence
                            
                            # Ensure proper orientation
                            if start_pos > end_pos:
                                start_pos, end_pos = end_pos, start_pos
                            
                            orf_length = end_pos - start_pos
                            
                            # Only include ORFs longer than 30 nucleotides (10 codons)
                            if orf_length >= 30:
                                orf_sequence = sequence[start_pos:end_pos]
                                orfs.append({
                                    'start': start_pos,
                                    'end': end_pos,
                                    'length': orf_length,
                                    'sequence': orf_sequence,
                                    'frame': -(frame + 1)  # Negative for reverse complement
                                })
                            break
        except Exception as e:
            print(f"DEBUG: Error processing reverse complement: {e}")
            # Continue without reverse complement analysis
        
        # Sort by length (largest first) and return the largest ORF
        if orfs:
            orfs.sort(key=lambda x: x['length'], reverse=True)
            largest_orf = orfs[0]
            
            print(f"DEBUG: Found {len(orfs)} ORFs")
            print(f"DEBUG: Largest ORF: {largest_orf['length']} bp, positions {largest_orf['start']}-{largest_orf['end']}, frame {largest_orf['frame']}")
            
            return largest_orf
        else:
            print("DEBUG: No ORFs found")
            return None
            
    except Exception as e:
        print(f"DEBUG: Error in find_largest_orf: {e}")
        return None

def find_3utr_region(sequence):
    """
    Find the 3' UTR region in the given DNA sequence.
    The 3' UTR is the region after the last stop codon of the largest ORF.
    Returns the start position, end position, and the 3' UTR sequence.
    """
    try:
        # First find the largest ORF to identify the end of coding sequence
        orf_info = find_largest_orf(sequence)
        
        if not orf_info:
            print("DEBUG: No ORF found, cannot identify 3' UTR")
            return None
        
        # The 3' UTR starts after the ORF ends
        utr_start = orf_info['end']
        utr_end = len(sequence)
        utr_length = utr_end - utr_start
        
        print(f"DEBUG: 3' UTR region: {utr_start}-{utr_end} ({utr_length} bp)")
        
        if utr_length < 20:
            print("DEBUG: 3' UTR too short, using full sequence")
            return None
        
        utr_sequence = sequence[utr_start:utr_end]
        
        return {
            'start': utr_start,
            'end': utr_end,
            'length': utr_length,
            'sequence': utr_sequence,
            'orf_end': orf_info['end']  # Keep track of where ORF ends
        }
        
    except Exception as e:
        print(f"DEBUG: Error in find_3utr_region: {e}")
        return None

def process_sequence_by_region(sequence, target_region, custom_start=None, custom_end=None):
    """Process sequence based on target region selection"""
    try:
        if target_region == 'full_gene':
            return sequence
        elif target_region == 'cds':
            # Find the largest ORF for CDS
            orf_info = find_largest_orf(sequence)
            if orf_info:
                print(f"DEBUG: Using CDS region: {orf_info['start']}-{orf_info['end']} ({orf_info['length']} bp)")
                # Return the full sequence with ORF info for highlighting
                return sequence, orf_info
            else:
                print("DEBUG: No ORF found, using full sequence")
                return sequence, None
        elif target_region == 'custom' and custom_start and custom_end:
            try:
                start = int(custom_start) - 1  # Convert to 0-based indexing
                end = int(custom_end)
                if 0 <= start < end <= len(sequence):
                    return sequence[start:end]
                else:
                    return None, "Custom coordinates out of range"
            except ValueError:
                return None, "Invalid custom coordinates"
        elif target_region == '3utr':
            # Find the 3' UTR region
            utr_info = find_3utr_region(sequence)
            if utr_info:
                print(f"DEBUG: Using 3' UTR region: {utr_info['start']}-{utr_info['end']} ({utr_info['length']} bp)")
                # Return the full sequence with UTR info for highlighting
                return sequence, utr_info
            else:
                print("DEBUG: No 3' UTR found, using full sequence")
                return sequence, None
        else:
            # For UTR, we'll use the full sequence for now
            # In a real implementation, you'd use Biopython to parse GenBank files
            # and extract the actual UTR regions
            return sequence
    except Exception as e:
        print(f"DEBUG: Error in process_sequence_by_region: {e}")
        # Fallback to full sequence
        return sequence, None

@app.route('/design', methods=['POST'])
def design_primers():
    """Design primers using primer3-py"""
    data = request.get_json()
    sequence = data.get('sequence', '').strip().upper()
    mode = data.get('mode', 'beginner')
    experiment_type = data.get('experiment_type', 'standard_pcr')
    target_region = data.get('target_region', 'full_gene')
    custom_start = data.get('custom_start')
    custom_end = data.get('custom_end')
    
    # For qPCR, always use 3' UTR as target region
    if experiment_type == 'qpcr':
        target_region = '3utr'
        print(f"DEBUG: qPCR experiment - automatically setting target region to 3' UTR")
    
    if not sequence:
        return jsonify({'error': 'Please enter a DNA sequence'})
    
    # Validate sequence
    if not re.match(r'^[ATCGN]+$', sequence):
        return jsonify({'error': 'Invalid DNA sequence. Please use only A, T, C, G, and N characters.'})
    
    # Process sequence based on target region
    processed_result = process_sequence_by_region(sequence, target_region, custom_start, custom_end)
    
    # Handle different return types
    if isinstance(processed_result, tuple):
        if len(processed_result) == 2:
            processed_sequence, orf_info = processed_result
        else:
            return jsonify({'error': processed_result[1]})
    else:
        processed_sequence = processed_result
        orf_info = None
    
    if isinstance(processed_sequence, tuple):
        return jsonify({'error': processed_sequence[1]})
    
    if len(processed_sequence) < 50:
        return jsonify({'error': 'Sequence too short. Please provide at least 50 nucleotides.'})
    
    try:
        # Configure primer3 parameters based on experiment type and target region
        if experiment_type == 'qpcr':
            # qPCR parameters: always target 3' UTR with small amplicons (<200 bp)
            if orf_info:
                # For 3' UTR qPCR, design primers to cover the 3' UTR region
                utr_length = orf_info['length']
                utr_start = orf_info['start']
                utr_end = orf_info['end']
                
                # qPCR: smaller amplicons for efficiency, cover 50-80% of UTR
                min_product_size = int(utr_length * 0.50)  # At least 50% of UTR
                max_product_size = int(utr_length * 0.80)  # Up to 80% of UTR
                
                # Ensure valid range and qPCR size limits (always <200bp)
                # First ensure min < max
                if min_product_size >= max_product_size:
                    min_product_size = max(80, utr_length // 3)
                    max_product_size = min(200, utr_length)
                
                # Then apply qPCR size constraints
                min_product_size = max(80, min_product_size)
                max_product_size = min(200, max_product_size)
                
                # Final validation: ensure min < max
                if min_product_size >= max_product_size:
                    # If still invalid, use safe defaults
                    min_product_size = 80
                    max_product_size = min(200, utr_length)
                    if min_product_size >= max_product_size:
                        max_product_size = min_product_size + 50  # Ensure at least 50bp difference
                
                print(f"DEBUG: qPCR 3' UTR - UTR: {utr_start}-{utr_end} ({utr_length}bp), Min: {min_product_size}, Max: {max_product_size}")
                
                # Use the 3' UTR sequence for primer design
                utr_sequence = sequence[utr_start:utr_end]
                
                primer3_config = {
                    'SEQUENCE_TEMPLATE': utr_sequence,
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 58.0,
                    'PRIMER_MAX_TM': 62.0,
                    'PRIMER_MIN_GC': 30.0,
                    'PRIMER_MAX_GC': 70.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[min_product_size, max_product_size]],
                    'PRIMER_MAX_HAIRPIN_TH': 2.0,
                    'PRIMER_MAX_SELF_ANY_TH': 2.0,
                    'PRIMER_MAX_SELF_END_TH': 2.0,
                    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 2.0,
                    'PRIMER_PAIR_MAX_COMPL_END_TH': 2.0
                }
            else:
                # Fallback: if no ORF found, use full sequence with qPCR constraints
                print(f"DEBUG: qPCR - No ORF found, using full sequence with qPCR constraints [80, 200]")
                primer3_config = {
                    'SEQUENCE_TEMPLATE': processed_sequence,
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 58.0,
                    'PRIMER_MAX_TM': 62.0,
                    'PRIMER_MIN_GC': 30.0,
                    'PRIMER_MAX_GC': 70.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[80, 200]],
                    'PRIMER_MAX_HAIRPIN_TH': 2.0,
                    'PRIMER_MAX_SELF_ANY_TH': 2.0,
                    'PRIMER_MAX_SELF_END_TH': 2.0,
                    'PRIMER_PAIR_MAX_COMPL_ANY_TH': 2.0,
                    'PRIMER_PAIR_MAX_COMPL_END_TH': 2.0
                }
        else:
            # Standard PCR parameters
            if target_region == 'full_gene':
                # For full gene, design primers to cover maximum of the gene
                sequence_length = len(sequence)
                
                # Standard PCR: Forward primer in first 50-100bp, reverse primer in last 50-100bp
                # Target 85-95% gene coverage
                min_product_size = int(sequence_length * 0.85)  # At least 85% of gene
                max_product_size = int(sequence_length * 0.95)  # Up to 95% of gene
                
                # Ensure valid range
                if min_product_size >= max_product_size:
                    min_product_size = int(sequence_length * 0.80)
                    max_product_size = int(sequence_length * 0.98)
                
                print(f"DEBUG: Standard PCR Full Gene - Sequence Length: {sequence_length}, Min: {min_product_size}, Max: {max_product_size}")
                
                # Use boundary constraints to position primers at beginning and end
                primer3_config = {
                    'SEQUENCE_TEMPLATE': sequence,
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[min_product_size, max_product_size]],
                    'PRIMER_LEFT_INPUT': [0, 100],  # Forward primer in first 100bp
                    'PRIMER_RIGHT_INPUT': [sequence_length-100, sequence_length]  # Reverse primer in last 100bp
                }
            elif target_region == 'cds' and orf_info:
                # For CDS standard PCR, ensure primers cover the entire ORF
                orf_length = orf_info['length']
                orf_start = orf_info['start']
                orf_end = orf_info['end']
                
                # Product size must be larger than ORF length and contain the entire ORF
                min_product_size = orf_length  # At least the same as ORF length
                max_product_size = orf_length + 200  # Allow extra nucleotides
                
                print(f"DEBUG: Standard PCR CDS - ORF: {orf_start}-{orf_end} ({orf_length}bp), Min: {min_product_size}, Max: {max_product_size}")
                
                # Forward primer must be before ORF start, reverse primer must be after ORF end
                # Use sequence masking to prevent primers in ORF region
                masked_sequence = sequence[:orf_start] + 'N' * (orf_end - orf_start) + sequence[orf_end:]
                print(f"DEBUG: Masked ORF region {orf_start}-{orf_end} for primer design")
                
                primer3_config = {
                    'SEQUENCE_TEMPLATE': masked_sequence,
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[min_product_size, max_product_size]]
                }
            elif target_region == '3utr' and orf_info:
                # For 3' UTR standard PCR, design primers to cover the 3' UTR region
                utr_length = orf_info['length']
                utr_start = orf_info['start']
                utr_end = orf_info['end']
                
                # Product size should cover most of the 3' UTR (80-95% coverage)
                min_product_size = int(utr_length * 0.70)  # At least 70% of UTR
                max_product_size = int(utr_length * 0.95)  # Up to 95% of UTR
                
                # Ensure valid range
                if min_product_size >= max_product_size:
                    min_product_size = max(200, utr_length // 2)
                    max_product_size = utr_length
                
                print(f"DEBUG: Standard PCR 3' UTR - UTR: {utr_start}-{utr_end} ({utr_length}bp), Min: {min_product_size}, Max: {max_product_size}")
                
                # Use the 3' UTR sequence for primer design
                utr_sequence = sequence[utr_start:utr_end]
                
                primer3_config = {
                    'SEQUENCE_TEMPLATE': utr_sequence,
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[min_product_size, max_product_size]]
                }
            elif target_region == 'custom' and custom_start and custom_end:
                # For custom coordinates, design primers within the specified region
                try:
                    start = int(custom_start) - 1  # Convert to 0-based indexing
                    end = int(custom_end)
                    custom_length = end - start
                    
                    if 0 <= start < end <= len(sequence):
                        # Product size should cover most of the custom region
                        min_product_size = int(custom_length * 0.80)  # At least 80% of custom region
                        max_product_size = custom_length  # Up to 100% of custom region
                        
                        print(f"DEBUG: Standard PCR Custom - Region: {start+1}-{end} ({custom_length}bp), Min: {min_product_size}, Max: {max_product_size}")
                        
                        # Use the custom region sequence for primer design
                        custom_sequence = sequence[start:end]
                        
                        primer3_config = {
                            'SEQUENCE_TEMPLATE': custom_sequence,
                            'PRIMER_TASK': 'generic',
                            'PRIMER_PICK_LEFT_PRIMER': 1,
                            'PRIMER_PICK_RIGHT_PRIMER': 1,
                            'PRIMER_OPT_SIZE': 20,
                            'PRIMER_MIN_SIZE': 18,
                            'PRIMER_MAX_SIZE': 25,
                            'PRIMER_OPT_TM': 60.0,
                            'PRIMER_MIN_TM': 57.0,
                            'PRIMER_MAX_TM': 63.0,
                            'PRIMER_MIN_GC': 20.0,
                            'PRIMER_MAX_GC': 80.0,
                            'PRIMER_PRODUCT_SIZE_RANGE': [[min_product_size, max_product_size]]
                        }
                    else:
                        return jsonify({'error': 'Custom coordinates out of range'})
                except ValueError:
                    return jsonify({'error': 'Invalid custom coordinates'})
            else:
                # Standard PCR parameters for other regions
                print(f"DEBUG: Standard PCR Other - Using default range [200, 1000]")
                primer3_config = {
                    'SEQUENCE_TEMPLATE': processed_sequence,
                    'PRIMER_TASK': 'generic',
                    'PRIMER_PICK_LEFT_PRIMER': 1,
                    'PRIMER_PICK_RIGHT_PRIMER': 1,
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_PRODUCT_SIZE_RANGE': [[200, 1000]]
                }
        
        # Run primer3
        result = primer3.bindings.design_primers(primer3_config, global_args={})
        
        # Check for errors
        if 'PRIMER_ERROR' in result and result['PRIMER_ERROR']:
            return jsonify({'error': f'Primer design failed: {result["PRIMER_ERROR"]}'})
        
        # Check if we got any primer pairs
        if 'PRIMER_PAIR' not in result or not result['PRIMER_PAIR']:
            return jsonify({'error': 'No suitable primer pairs found. Try a longer sequence or different parameters.'})
        
        # Extract primer pairs with detailed analysis
        primer_pairs = []
        num_pairs = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
        
        for i in range(min(3, num_pairs)):
            # Get primer sequences from individual keys
            forward_primer = result.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')
            reverse_primer = result.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')
            product_size = result.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0)
            
            if not forward_primer or not reverse_primer:
                continue
            
            # Calculate detailed analysis
            score, explanations = calculate_primer_score(forward_primer, reverse_primer, product_size, experiment_type)
            
            # Get detailed primer properties
            forward_tm = round(primer3.calc_tm(forward_primer), 1)
            reverse_tm = round(primer3.calc_tm(reverse_primer), 1)
            forward_gc = round((forward_primer.count('G') + forward_primer.count('C')) / len(forward_primer) * 100, 1)
            reverse_gc = round((reverse_primer.count('G') + reverse_primer.count('C')) / len(reverse_primer) * 100, 1)
            
            # Hairpin analysis
            forward_hairpin = primer3.calc_hairpin(forward_primer)
            reverse_hairpin = primer3.calc_hairpin(reverse_primer)
            
            # Self-dimer analysis
            forward_dimer = primer3.calc_homodimer(forward_primer)
            reverse_dimer = primer3.calc_homodimer(reverse_primer)
            
            # Hetero-dimer analysis
            hetero_dimer = primer3.calc_heterodimer(forward_primer, reverse_primer)
            
            # GC clamp analysis
            forward_gc_clamp = forward_primer[-1] in 'GC'
            reverse_gc_clamp = reverse_primer[-1] in 'GC'
            
            # Analyze secondary structures
            forward_secondary = analyze_secondary_structures(forward_primer)
            reverse_secondary = analyze_secondary_structures(reverse_primer)
            
            # Generate dimer visualizations
            forward_self_dimer = visualize_dimer(forward_primer, forward_primer, "self")
            reverse_self_dimer = visualize_dimer(reverse_primer, reverse_primer, "self")
            hetero_dimer_viz = visualize_dimer(forward_primer, reverse_primer, "hetero")
            
            primer_pairs.append({
                'forward': forward_primer,
                'reverse': reverse_primer,
                'product_size': product_size,
                'score': score,
                'score_label': get_score_label(score),
                'explanations': explanations,
                'forward_tm': forward_tm,
                'reverse_tm': reverse_tm,
                'forward_gc': forward_gc,
                'reverse_gc': reverse_gc,
                'forward_length': len(forward_primer),
                'reverse_length': len(reverse_primer),
                'forward_gc_clamp': forward_gc_clamp,
                'reverse_gc_clamp': reverse_gc_clamp,
                'forward_hairpin_dg': round(forward_hairpin.dg, 2),
                'reverse_hairpin_dg': round(reverse_hairpin.dg, 2),
                'forward_dimer_dg': round(forward_dimer.dg, 2),
                'reverse_dimer_dg': round(reverse_dimer.dg, 2),
                'hetero_dimer_dg': round(hetero_dimer.dg, 2),
                'hairpin_risky': forward_hairpin.dg < -2 or reverse_hairpin.dg < -2,
                'dimer_risky': forward_dimer.dg < -6 or reverse_dimer.dg < -6 or hetero_dimer.dg < -7,
                'summary': get_primer_summary(score, explanations),
                'forward_secondary': forward_secondary,
                'reverse_secondary': reverse_secondary,
                'forward_self_dimer': forward_self_dimer,
                'reverse_self_dimer': reverse_self_dimer,
                'hetero_dimer_viz': hetero_dimer_viz
            })
        
        # Sort by score (highest first) and take top 3
        primer_pairs.sort(key=lambda x: x['score'], reverse=True)
        top_primers = primer_pairs[:3]
        
        # Add primer positions for sequence view and apply full gene bonus
        for i, primer in enumerate(top_primers):
            try:
                # Get primer positions from primer3 output - handle both list and int formats
                forward_pos_data = result.get(f'PRIMER_LEFT_{i}', 0)
                reverse_pos_data = result.get(f'PRIMER_RIGHT_{i}', 0)
                
                # Extract position values safely
                if isinstance(forward_pos_data, list) and len(forward_pos_data) > 0:
                    forward_pos = forward_pos_data[0]
                else:
                    forward_pos = forward_pos_data
                    
                if isinstance(reverse_pos_data, list) and len(reverse_pos_data) > 0:
                    reverse_pos = reverse_pos_data[0]
                else:
                    reverse_pos = reverse_pos_data
                
                # Ensure positions are integers
                forward_pos = int(forward_pos) if forward_pos is not None else 0
                reverse_pos = int(reverse_pos) if reverse_pos is not None else 0
                
                # Forward primer positions (direct sequence match)
                if target_region == 'cds' and orf_info:
                    # For CDS, primers are designed on full sequence with constraints
                    # No need for ORF offset adjustment
                    primer['forward_start'] = forward_pos
                    primer['forward_end'] = forward_pos + len(primer['forward'])
                else:
                    primer['forward_start'] = forward_pos
                    primer['forward_end'] = forward_pos + len(primer['forward'])
                
                # Reverse primer positions (binds to reverse complement)
                # The reverse primer binds to the reverse complement, so we need to calculate
                # where the reverse complement of the reverse primer appears in the sequence
                reverse_complement = get_complement(primer['reverse'])[::-1]  # Reverse complement
                
                # Find where the reverse complement appears in the sequence
                # This is where the reverse primer actually binds
                if target_region == 'cds' and orf_info:
                    # For CDS, search in the full sequence
                    reverse_binding_pos = sequence.find(reverse_complement)
                    if reverse_binding_pos != -1:
                        primer['reverse_start'] = reverse_binding_pos
                        primer['reverse_end'] = reverse_binding_pos + len(primer['reverse'])
                    else:
                        # Fallback: use the original reverse_pos
                        primer['reverse_start'] = reverse_pos
                        primer['reverse_end'] = reverse_pos + len(primer['reverse'])
                else:
                    reverse_binding_pos = processed_sequence.find(reverse_complement)
                    if reverse_binding_pos != -1:
                        primer['reverse_start'] = reverse_binding_pos
                        primer['reverse_end'] = reverse_binding_pos + len(primer['reverse'])
                    else:
                        # Fallback: use the original reverse_pos but adjust for reverse complement
                        primer['reverse_start'] = reverse_pos
                        primer['reverse_end'] = reverse_pos + len(primer['reverse'])
                
                # Apply full gene coverage bonus for "Full Gene" target region
                if target_region == 'full_gene':
                    # Calculate coverage percentage
                    total_sequence_length = len(processed_sequence)
                    coverage_length = primer['reverse_end'] - primer['forward_start']
                    coverage_percentage = (coverage_length / total_sequence_length) * 100
                    
                    # Enhanced scoring for Standard PCR - prioritize 85-95% coverage
                    if experiment_type == 'standard_pcr':
                        if coverage_percentage >= 95:
                            primer['score'] = min(100, primer['score'] + 20)  # Maximum bonus
                            primer['full_gene_coverage'] = f"Excellent coverage ({coverage_percentage:.1f}% of gene)"
                        elif coverage_percentage >= 90:
                            primer['score'] = min(100, primer['score'] + 15)
                            primer['full_gene_coverage'] = f"Very good coverage ({coverage_percentage:.1f}% of gene)"
                        elif coverage_percentage >= 85:
                            primer['score'] = min(100, primer['score'] + 10)
                            primer['full_gene_coverage'] = f"Good coverage ({coverage_percentage:.1f}% of gene)"
                        elif coverage_percentage >= 75:
                            primer['score'] = min(100, primer['score'] + 5)
                            primer['full_gene_coverage'] = f"Moderate coverage ({coverage_percentage:.1f}% of gene)"
                        else:
                            primer['full_gene_coverage'] = f"Limited coverage ({coverage_percentage:.1f}% of gene) - Consider alternatives"
                    else:
                        # Original scoring for qPCR
                        if coverage_percentage >= 80:
                            primer['score'] = min(100, primer['score'] + 15)  # Maximum bonus
                            primer['full_gene_coverage'] = f"Excellent coverage ({coverage_percentage:.1f}% of gene)"
                        elif coverage_percentage >= 60:
                            primer['score'] = min(100, primer['score'] + 10)
                            primer['full_gene_coverage'] = f"Good coverage ({coverage_percentage:.1f}% of gene)"
                        elif coverage_percentage >= 40:
                            primer['score'] = min(100, primer['score'] + 5)
                            primer['full_gene_coverage'] = f"Moderate coverage ({coverage_percentage:.1f}% of gene)"
                        else:
                            primer['full_gene_coverage'] = f"Limited coverage ({coverage_percentage:.1f}% of gene)"
                    
                    # Update score label after bonus
                    primer['score_label'] = get_score_label(primer['score'])
                    
                    print(f"DEBUG: Full gene coverage for primer {i+1}: {coverage_percentage:.1f}%")
                
                # Calculate CDS coverage for CDS target region
                if target_region == 'cds' and orf_info:
                    cds_coverage = calculate_cds_coverage(
                        primer['forward_start'], 
                        primer['reverse_end'], 
                        orf_info
                    )
                    primer['cds_coverage'] = f"{cds_coverage}% CDS coverage"
                    print(f"DEBUG: CDS coverage for primer {i+1}: {cds_coverage}%")
                    
                    # Add bonus points for better CDS coverage (90-100% target)
                    if cds_coverage >= 100:
                        primer['score'] = min(100, primer['score'] + 15)  # Maximum bonus
                        primer['cds_coverage_label'] = f"Perfect CDS coverage ({cds_coverage}%)"
                    elif cds_coverage >= 95:
                        primer['score'] = min(100, primer['score'] + 12)
                        primer['cds_coverage_label'] = f"Excellent CDS coverage ({cds_coverage}%)"
                    elif cds_coverage >= 90:
                        primer['score'] = min(100, primer['score'] + 10)
                        primer['cds_coverage_label'] = f"Good CDS coverage ({cds_coverage}%)"
                    elif cds_coverage >= 80:
                        primer['score'] = min(100, primer['score'] + 5)
                        primer['cds_coverage_label'] = f"Moderate CDS coverage ({cds_coverage}%)"
                    else:
                        primer['cds_coverage_label'] = f"Limited CDS coverage ({cds_coverage}%) - Consider alternatives"
                    
                    # Update score label after bonus
                    primer['score_label'] = get_score_label(primer['score'])
                
                # Calculate 3' UTR coverage for 3' UTR target region
                if target_region == '3utr' and orf_info:
                    # For 3' UTR, orf_info actually contains UTR information
                    # Adjust positions for 3' UTR (primer positions are relative to UTR sequence)
                    utr_start = orf_info['start']
                    utr_end = orf_info['end']
                    
                    # For 3' UTR, we need to handle reverse primer binding differently
                    # The reverse primer binds to the reverse complement of the UTR sequence
                    utr_sequence = sequence[utr_start:utr_end]
                    reverse_complement = get_complement(primer['reverse'])[::-1]
                    
                    # Find where the reverse complement appears in the UTR sequence
                    reverse_binding_pos_in_utr = utr_sequence.find(reverse_complement)
                    
                    if reverse_binding_pos_in_utr != -1:
                        # Convert UTR-relative positions to full sequence positions
                        full_forward_start = utr_start + primer['forward_start']
                        full_forward_end = utr_start + primer['forward_end']
                        full_reverse_start = utr_start + reverse_binding_pos_in_utr
                        full_reverse_end = utr_start + reverse_binding_pos_in_utr + len(primer['reverse'])
                        
                        # Update primer positions to full sequence coordinates
                        primer['forward_start'] = full_forward_start
                        primer['forward_end'] = full_forward_end
                        primer['reverse_start'] = full_reverse_start
                        primer['reverse_end'] = full_reverse_end
                    else:
                        # Fallback: use the original positions but adjust for UTR offset
                        full_forward_start = utr_start + primer['forward_start']
                        full_forward_end = utr_start + primer['forward_end']
                        full_reverse_start = utr_start + primer['reverse_start']
                        full_reverse_end = utr_start + primer['reverse_end']
                        
                        # Update primer positions to full sequence coordinates
                        primer['forward_start'] = full_forward_start
                        primer['forward_end'] = full_forward_end
                        primer['reverse_start'] = full_reverse_start
                        primer['reverse_end'] = full_reverse_end
                    
                    # Calculate UTR coverage
                    utr_length = orf_info['length']
                    coverage_length = primer['reverse_end'] - primer['forward_start']
                    utr_coverage = (coverage_length / utr_length) * 100
                    
                    primer['utr_coverage'] = f"{utr_coverage:.1f}% 3' UTR coverage"
                    print(f"DEBUG: 3' UTR coverage for primer {i+1}: {utr_coverage:.1f}%")
                    print(f"DEBUG: 3' UTR reverse binding - UTR sequence: {utr_sequence[:50]}...")
                    print(f"DEBUG: 3' UTR reverse binding - Reverse complement: {reverse_complement}")
                    print(f"DEBUG: 3' UTR reverse binding - Position in UTR: {reverse_binding_pos_in_utr}")
                    
                    # Add bonus points for better UTR coverage
                    if utr_coverage >= 90:
                        primer['score'] = min(100, primer['score'] + 10)  # Maximum bonus
                        primer['utr_coverage_label'] = f"Excellent 3' UTR coverage ({utr_coverage:.1f}%)"
                    elif utr_coverage >= 80:
                        primer['score'] = min(100, primer['score'] + 7)
                        primer['utr_coverage_label'] = f"Good 3' UTR coverage ({utr_coverage:.1f}%)"
                    elif utr_coverage >= 70:
                        primer['score'] = min(100, primer['score'] + 5)
                        primer['utr_coverage_label'] = f"Moderate 3' UTR coverage ({utr_coverage:.1f}%)"
                    else:
                        primer['utr_coverage_label'] = f"Limited 3' UTR coverage ({utr_coverage:.1f}%)"
                    
                    # Update score label after bonus
                    primer['score_label'] = get_score_label(primer['score'])
                
                print(f"DEBUG: Primer {i+1} positions:")
                print(f"  Forward: {primer['forward_start']}-{primer['forward_end']} ({primer['forward']})")
                print(f"  Reverse: {primer['reverse_start']}-{primer['reverse_end']} ({primer['reverse']})")
                print(f"  Reverse complement: {reverse_complement}")
                print(f"  Reverse binding position: {reverse_binding_pos}")
                
            except Exception as e:
                print(f"ERROR calculating primer positions: {e}")
                # Fallback to simple positioning if there's an error
                primer['forward_start'] = 100 + (i * 50)
                primer['forward_end'] = 100 + (i * 50) + len(primer['forward'])
                primer['reverse_start'] = 500 + (i * 50)
                primer['reverse_end'] = 500 + (i * 50) + len(primer['reverse'])
        
        # Re-sort by updated scores for full gene coverage
        if target_region == 'full_gene':
            top_primers.sort(key=lambda x: x['score'], reverse=True)
        
        # Create individual highlighted sequences for each primer set
        individual_highlighted_sequences = []
        for i, primer in enumerate(top_primers):
            try:
                individual_sequence = get_highlighted_sequence(processed_sequence, [primer], i, orf_info, None)
                individual_highlighted_sequences.append(individual_sequence)
            except Exception as e:
                print(f"DEBUG: Error creating highlighted sequence for primer {i+1}: {e}")
                fallback_sequence = [{'char': char, 'color_class': ''} for char in processed_sequence]
                individual_highlighted_sequences.append(fallback_sequence)
        
        # Create the main highlighted sequence (all primers)
        try:
            highlighted_sequence = get_highlighted_sequence(processed_sequence, top_primers, 0, orf_info, None)
        except Exception as e:
            print(f"DEBUG: Error creating main highlighted sequence: {e}")
            highlighted_sequence = [{'char': char, 'color_class': ''} for char in processed_sequence]
        
        # Generate unique result ID
        result_id = str(uuid.uuid4())
        
        # Store results in memory
        result_store[result_id] = {
            'primers': top_primers,
            'sequence': processed_sequence,
            'original_sequence': sequence,
            'sequence_length': len(processed_sequence),
            'experiment_type': experiment_type,
            'mode': mode,
            'orf_info': orf_info,
            'utr_info': None,
            'highlighted_sequence': highlighted_sequence,
            'individual_highlighted_sequences': individual_highlighted_sequences,
            'target_region': target_region,
            'custom_start': custom_start,
            'custom_end': custom_end
        }
        
        # Return redirect URL with short result ID
        return jsonify({
            'redirect': f'/result/{result_id}'
        })
        
    except Exception as e:
        return jsonify({'error': f'Error designing primers: {str(e)}'})

def get_primer_summary(score, explanations):
    """Generate a one-line summary of primer quality"""
    if score >= 80:
        return "Excellent primer pair with optimal parameters"
    elif score >= 60:
        return "Good primer pair with minor issues"
    else:
        return "Primer pair has significant issues - consider alternatives"

@app.route('/download-csv', methods=['POST'])
def download_csv():
    """Download primers as CSV"""
    data = request.get_json()
    primers = data.get('primers', [])
    
    if not primers:
        return jsonify({'error': 'No primers to download'})
    
    # Create CSV in memory
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Write header
    writer.writerow(['Primer Pair', 'Forward Primer', 'Reverse Primer', 'Product Size (bp)', 
                    'Score', 'Forward Tm (°C)', 'Reverse Tm (°C)', 'Forward GC (%)', 'Reverse GC (%)'])
    
    # Write data
    for i, primer in enumerate(primers, 1):
        writer.writerow([
            f'Pair {i}',
            primer['forward'],
            primer['reverse'],
            primer['product_size'],
            f"{primer['score']}/100",
            primer['forward_tm'],
            primer['reverse_tm'],
            primer['forward_gc'],
            primer['reverse_gc']
        ])
    
    output.seek(0)
    
    # Create response
    response = app.response_class(
        output.getvalue(),
        mimetype='text/csv',
        headers={'Content-Disposition': f'attachment; filename=primers_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'}
    )
    
    return response

def analyze_secondary_structures(primer):
    """Analyze secondary structures and return visualization data"""
    # Calculate hairpin potential
    hairpin = primer3.calc_hairpin(primer)
    
    # Calculate self-dimer potential
    dimer = primer3.calc_homodimer(primer)
    
    # Simple hairpin detection (look for palindromic regions)
    hairpin_region = None
    if len(primer) >= 6:
        for i in range(len(primer) - 5):
            for j in range(i + 6, len(primer) - 1):
                region1 = primer[i:i+3]
                region2 = primer[j:j+3]
                # Check if regions are complementary
                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                if region2 == ''.join(complement.get(base, base) for base in region1[::-1]):
                    hairpin_region = [i, j, region1, region2]  # Use list instead of tuple
                    break
            if hairpin_region:
                break
    
    return {
        'hairpin_dg': hairpin.dg,
        'dimer_dg': dimer.dg,
        'hairpin_risky': hairpin.dg < -2,
        'dimer_risky': dimer.dg < -6,
        'hairpin_region': hairpin_region
    }

def get_complement(sequence):
    """Get complementary sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in sequence)

@app.route('/test')
def test():
    """Test route to check if template rendering works"""
    try:
        return render_template('results.html', 
                             primers=[], 
                             sequence_length=0,
                             experiment_type='standard_pcr',
                             sequence='',
                             highlighted_sequence=[],
                             individual_highlighted_sequences=[],
                             orf_info=None,
                             utr_info=None)
    except Exception as e:
        import traceback
        traceback.print_exc()
        return f"Template error: {str(e)}", 500

def calculate_cds_coverage(forward_pos, reverse_pos, orf_info):
    """Calculate what percentage of the CDS is covered by the primer pair"""
    if not orf_info:
        return 100.0  # Default if no ORF info
    
    orf_start = orf_info['start']
    orf_end = orf_info['end']
    orf_length = orf_info['length']
    
    # Calculate the amplicon region
    amplicon_start = forward_pos
    amplicon_end = reverse_pos + 20  # Add primer length for end position
    
    # Calculate overlap with ORF
    overlap_start = max(amplicon_start, orf_start)
    overlap_end = min(amplicon_end, orf_end)
    
    if overlap_end <= overlap_start:
        return 0.0  # No overlap
    
    overlap_length = overlap_end - overlap_start
    coverage_percentage = (overlap_length / orf_length) * 100
    
    return round(coverage_percentage, 1)

@app.route('/version')
def version():
    """Return version information"""
    return jsonify({
        'version': __version__,
        'author': __author__,
        'organization': __organization__,
        'license': __license__,
        'copyright': __copyright__,
        'website': __website__,
        'description': 'PCR Primer Design Tool'
    })

@app.route('/learn')
def learn():
    """Display primer design theory and educational content"""
    return render_template('learn.html')

import os

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))  # 5000 is fallback for local
    app.run(host='0.0.0.0', port=port, debug=True)