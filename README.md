# Primerly - PCR Primer Design Tool

**A Biovagon Initiative**

A beginner-friendly web application for designing PCR primers using primer3-py and Biopython. Primerly provides an intuitive interface for both standard PCR and qPCR primer design with comprehensive quality analysis.

## 🌟 Features

### Core Functionality
- **Dual Experiment Types**: Standard PCR and qPCR with optimized parameters
- **Multiple Input Methods**: Direct sequence input or NCBI accession number fetching
- **Target Region Selection**: Full gene, CDS (with ORF detection), 3' UTR, or custom coordinates
- **Comprehensive Analysis**: Tm, GC content, hairpins, dimers, secondary structures
- **Visual Feedback**: Color-coded quality indicators and interactive gene maps
- **Export Options**: CSV download with detailed primer information

### User Experience
- **Beginner/Expert Modes**: Simplified interface with progressive disclosure
- **Mobile-Friendly Design**: Responsive Bootstrap-based interface
- **Educational Content**: Built-in primer design theory and explanations
- **Real-time Validation**: Immediate feedback on input parameters

## 🚀 Quick Start

### Prerequisites
- Python 3.8 or higher
- pip package manager

### Installation
```bash
# Clone the repository
git clone https://github.com/biovagon/primerly.git
cd primerly

# Install dependencies
pip install -r requirements.txt

# Run the application
python app.py
```

### Usage
1. Open your browser and navigate to `https://primerly.onrender.com/`
2. Choose your experiment type (Standard PCR or qPCR)
3. Select your mode (Beginner or Expert)
4. Enter your sequence or NCBI accession number
5. Choose your target region
6. Click "Design Primers" to get results

## 📋 Requirements

- Flask==2.3.3
- primer3-py==0.6.1
- biopython==1.81
- gunicorn==21.2.0 (for production)

## 🏗️ Architecture

### Backend
- **Flask**: Web framework for handling HTTP requests
- **primer3-py**: Core primer design algorithm
- **Biopython**: Sequence analysis and NCBI integration
- **CSV Export**: Results download functionality

### Frontend
- **Bootstrap 5**: Responsive UI framework
- **JavaScript**: Interactive features and validation
- **HTML5/CSS3**: Modern web standards

## 📁 Project Structure

```
Primerly/
├── app.py                 # Main Flask application
├── requirements.txt       # Python dependencies
├── templates/            # HTML templates
│   ├── index.html        # Homepage
│   ├── results.html      # Results page
│   └── learn.html        # Educational content
├── Procfile             # Heroku deployment
├── wsgi.py              # WSGI entry point
├── desktop_app.py       # Desktop application wrapper
├── build_app.py         # PyInstaller build script
└── DEPLOYMENT_GUIDE.md  # Deployment instructions
```

## 🎯 Target Users

- **Students**: Learning PCR primer design concepts
- **Researchers**: Quick primer design for experiments
- **Educators**: Teaching molecular biology concepts
- **Biotech Professionals**: Rapid primer validation

## 🔬 Scientific Background

Primerly implements industry-standard primer design algorithms based on:
- **Melting Temperature (Tm)**: Optimal range 55-65°C
- **GC Content**: Target 40-60% for stability
- **GC Clamp**: 3' end stability for efficient extension
- **Secondary Structures**: Hairpin and dimer analysis
- **Product Size**: Optimized for different experiment types

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup
```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install development dependencies
pip install -r requirements.txt

# Run in development mode
python app.py
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👨‍💻 Author

**Jeevitha C M** - Founder of Biovagon
- Website: [https://www.biovagon.org/](https://www.biovagon.org/)
- Email: jeevithacm21@gmail.com

## 🙏 Acknowledgments

- **primer3-py**: Core primer design algorithms
- **Biopython**: Sequence analysis tools
- **Flask**: Web framework
- **Bootstrap**: UI framework
- **Scientific Community**: For feedback and testing

## 📞 Support

For support, questions, or feature requests:
- Create an issue on GitHub
- Contact: jeevithacm21@gmail.com
- Visit: [https://www.biovagon.org/](https://www.biovagon.org/)

---

**Primerly** - Empowering Life Sciences Enthusiasts 🌱💡
*A Biovagon Initiative* 
