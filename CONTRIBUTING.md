# Contributing to Primerly

Thank you for your interest in contributing to Primerly! This document provides guidelines and information for contributors.

## 🤝 How to Contribute

### Types of Contributions

We welcome various types of contributions:

- **Bug Reports**: Report issues and bugs
- **Feature Requests**: Suggest new features and improvements
- **Code Contributions**: Submit pull requests with code changes
- **Documentation**: Improve documentation and tutorials
- **Testing**: Help test the application and report issues
- **Translation**: Help translate the interface to other languages

## 🐛 Reporting Bugs

### Before Submitting a Bug Report

1. **Check existing issues**: Search the issue tracker to see if the bug has already been reported
2. **Test with latest version**: Ensure you're using the latest version of Primerly
3. **Reproduce the issue**: Try to reproduce the bug consistently

### Bug Report Template

When reporting a bug, please include:

```markdown
**Bug Description**
A clear and concise description of the bug.

**Steps to Reproduce**
1. Go to '...'
2. Click on '...'
3. Enter '...'
4. See error

**Expected Behavior**
What you expected to happen.

**Actual Behavior**
What actually happened.

**Environment**
- OS: [e.g. Windows 10, macOS, Linux]
- Python Version: [e.g. 3.8, 3.9, 3.10]
- Browser: [e.g. Chrome, Firefox, Safari]
- Primerly Version: [e.g. 1.0.0]

**Additional Information**
Any other context, screenshots, or error messages.
```

## 💡 Suggesting Features

### Feature Request Template

```markdown
**Feature Description**
A clear and concise description of the feature.

**Use Case**
Why this feature would be useful and how it would be used.

**Proposed Implementation**
Any ideas on how this could be implemented (optional).

**Alternatives Considered**
Any alternative solutions you've considered (optional).
```

## 🔧 Code Contributions

### Development Setup

1. **Fork the repository**
   ```bash
   git clone https://github.com/your-username/Primerly.git
   cd Primerly
   ```

2. **Create a virtual environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

4. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

### Coding Standards

- **Python**: Follow PEP 8 style guidelines
- **JavaScript**: Use consistent indentation and naming conventions
- **HTML/CSS**: Follow Bootstrap conventions
- **Comments**: Add clear, descriptive comments for complex logic
- **Documentation**: Update relevant documentation for new features

### Testing

- **Test your changes**: Ensure the application works correctly
- **Test edge cases**: Consider unusual inputs and scenarios
- **Cross-browser testing**: Test on different browsers if UI changes
- **Mobile testing**: Verify mobile responsiveness

### Pull Request Process

1. **Update documentation**: Update README, CHANGELOG, or other docs as needed
2. **Add tests**: Include tests for new functionality if applicable
3. **Update requirements**: Add new dependencies to requirements.txt if needed
4. **Submit PR**: Create a pull request with a clear description

### Pull Request Template

```markdown
**Description**
Brief description of changes made.

**Type of Change**
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

**Testing**
- [ ] Tested locally
- [ ] Added/updated tests
- [ ] All tests pass

**Screenshots**
If applicable, add screenshots to help explain your changes.

**Checklist**
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No new warnings generated
```

## 📚 Documentation Contributions

### Areas for Documentation

- **User Guides**: Tutorials and how-to guides
- **API Documentation**: Code documentation and examples
- **Installation Guides**: Setup instructions for different platforms
- **Troubleshooting**: Common issues and solutions

### Documentation Standards

- Use clear, concise language
- Include code examples where appropriate
- Add screenshots for UI-related documentation
- Keep documentation up-to-date with code changes

## 🧪 Testing Contributions

### Testing Areas

- **Unit Tests**: Test individual functions and components
- **Integration Tests**: Test how components work together
- **User Interface Tests**: Test the web interface
- **Performance Tests**: Test application performance
- **Cross-platform Tests**: Test on different operating systems

### Running Tests

```bash
# Run all tests
python -m pytest

# Run specific test file
python -m pytest test_app.py

# Run with coverage
python -m pytest --cov=app
```

## 🌍 Translation Contributions

### Supported Languages

We're working on supporting multiple languages. Currently:
- English (primary)
- Spanish (in progress)
- French (planned)
- German (planned)

### Translation Process

1. **Identify strings**: Find translatable text in templates and code
2. **Create translation files**: Add translations to language files
3. **Test translations**: Verify translations work correctly
4. **Submit changes**: Create pull request with translation updates

## 📋 Code of Conduct

### Our Standards

- **Be respectful**: Treat all contributors with respect
- **Be inclusive**: Welcome contributors from diverse backgrounds
- **Be constructive**: Provide constructive feedback
- **Be patient**: Understand that contributors have different skill levels

### Unacceptable Behavior

- Harassment or discrimination
- Trolling or insulting comments
- Publishing others' private information
- Other conduct inappropriate in a professional setting

## 🏆 Recognition

### Contributors

We recognize contributors in several ways:

- **Contributors list**: All contributors are listed in the README
- **Release notes**: Contributors are credited in release notes
- **Special thanks**: Significant contributors receive special recognition

### Contribution Levels

- **Bronze**: 1-5 contributions
- **Silver**: 6-15 contributions
- **Gold**: 16+ contributions
- **Platinum**: Major contributions or long-term involvement

## 📞 Getting Help

### Communication Channels

- **GitHub Issues**: For bug reports and feature requests
- **GitHub Discussions**: For general questions and discussions
- **Email**: For private or sensitive matters

### Resources

- **Documentation**: Check the README and other docs first
- **Code Examples**: Look at existing code for patterns
- **Community**: Ask questions in GitHub Discussions

## 🎯 Project Goals

### Current Focus Areas

- **Performance**: Improve application speed and efficiency
- **Usability**: Enhance user experience and interface
- **Reliability**: Increase stability and error handling
- **Documentation**: Improve and expand documentation
- **Testing**: Increase test coverage and quality

### Long-term Vision

- **Accessibility**: Make the application accessible to all users
- **Internationalization**: Support for multiple languages
- **Mobile Apps**: Native mobile applications
- **Advanced Features**: Machine learning and AI integration

---

Thank you for contributing to Primerly! Your contributions help make PCR primer design accessible to everyone. 🧬 