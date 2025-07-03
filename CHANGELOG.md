# Changelog

All notable changes to Primerly will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-12-29

### Added
- **Core Application**: Complete Flask web application for PCR primer design
- **Dual Experiment Types**: Support for Standard PCR and qPCR with optimized parameters
- **Multiple Input Methods**: Direct sequence input and NCBI accession number fetching
- **Target Region Selection**: Full gene, CDS (with ORF detection), 3' UTR, and custom coordinates
- **Comprehensive Quality Analysis**: Tm, GC content, hairpins, dimers, secondary structures
- **Visual Feedback**: Color-coded quality indicators and interactive gene maps
- **Beginner/Expert Modes**: Simplified interface with progressive disclosure
- **Mobile-Friendly Design**: Responsive Bootstrap-based interface
- **Educational Content**: Built-in primer design theory and explanations
- **Export Functionality**: CSV download with detailed primer information
- **Deployment Options**: Heroku, PythonAnywhere, and desktop application support
- **Copyright and Branding**: Updated to reflect Biovagon initiative

### Technical Features
- **ORF Detection**: Automatic identification of largest open reading frame
- **Primer Positioning Logic**: Intelligent placement based on target region
- **Quality Scoring**: Comprehensive evaluation with visual indicators
- **Secondary Structure Analysis**: Hairpin and dimer visualization
- **Gene Map Visualization**: Interactive display of primer binding positions
- **Error Handling**: Robust validation and user-friendly error messages

### User Experience
- **Intuitive Interface**: Clean, modern design with clear navigation
- **Real-time Validation**: Immediate feedback on input parameters
- **Educational Resources**: Comprehensive learning materials
- **Export Options**: Multiple formats for result sharing

### Documentation
- **README.md**: Comprehensive project documentation
- **DEPLOYMENT_GUIDE.md**: Step-by-step deployment instructions
- **LICENSE**: MIT License for open source distribution
- **CHANGELOG.md**: Version history and updates tracking

### Deployment
- **Web Deployment**: Heroku and PythonAnywhere configuration files
- **Desktop Application**: PyInstaller and Tkinter-based desktop app
- **Requirements**: Complete dependency management

## [0.9.0] - 2024-12-28

### Added
- Initial development version
- Basic primer design functionality
- Simple web interface
- Core primer3-py integration

### Changed
- Multiple iterations of UI/UX improvements
- Enhanced primer quality analysis
- Improved error handling
- Performance optimizations

### Fixed
- Various bug fixes and stability improvements
- Syntax errors and compatibility issues
- Template rendering problems
- Data validation issues

---

## Version History Summary

### Major Features by Version
- **v1.0.0**: Production-ready application with all core features
- **v0.9.0**: Development and testing phase

### Key Milestones
- ✅ **Core Functionality**: Complete primer design pipeline
- ✅ **User Interface**: Professional, responsive web design
- ✅ **Quality Analysis**: Comprehensive scoring and visualization
- ✅ **Educational Content**: Built-in learning resources
- ✅ **Deployment Ready**: Multiple deployment options
- ✅ **Documentation**: Complete project documentation

### Future Roadmap
- **v1.1.0**: Batch processing capabilities
- **v1.2.0**: Advanced algorithms and machine learning
- **v1.3.0**: Cloud integration and API access
- **v2.0.0**: Mobile applications and enhanced features

---

For detailed information about each version, please refer to the commit history and release notes.

---

**Author**: Jeevitha C M, Founder of Biovagon  
**Organization**: Biovagon  
**Website**: https://www.biovagon.org/  
**License**: MIT 