# Primerly Deployment Guide

**A Biovagon Initiative**

This guide provides comprehensive instructions for deploying Primerly, a PCR primer design tool developed by Jeevitha C M, Founder of Biovagon.

## Option 1: Deploy as Web Application

### A. Deploy to Heroku (Recommended for beginners)

**Prerequisites:**
- Heroku account (free at heroku.com)
- Git installed
- Heroku CLI installed

**Steps:**
1. **Install Heroku CLI:**
   ```bash
   # Windows
   winget install --id=Heroku.HerokuCLI
   
   # Mac
   brew tap heroku/brew && brew install heroku
   ```

2. **Login to Heroku:**
   ```bash
   heroku login
   ```

3. **Create Heroku app:**
   ```bash
   heroku create your-primerly-app-name
   ```

4. **Deploy:**
   ```bash
   git add .
   git commit -m "Initial deployment"
   git push heroku main
   ```

5. **Open your app:**
   ```bash
   heroku open
   ```

### B. Deploy to PythonAnywhere (Free tier available)

**Steps:**
1. Create account at pythonanywhere.com
2. Upload your files via Files tab
3. Create a new web app
4. Set the WSGI file to point to `wsgi.py`
5. Configure your domain

### C. Deploy to VPS (Advanced)

**Steps:**
1. **Set up server:**
   ```bash
   sudo apt update
   sudo apt install python3-pip nginx
   ```

2. **Install dependencies:**
   ```bash
   pip3 install -r requirements.txt
   ```

3. **Configure Nginx:**
   ```bash
   sudo nano /etc/nginx/sites-available/primerly
   ```

4. **Set up systemd service:**
   ```bash
   sudo nano /etc/systemd/system/primerly.service
   ```

## Option 2: Convert to Desktop Application

### A. Using PyInstaller (Windows/Mac/Linux)

**Prerequisites:**
- PyInstaller installed: `pip install pyinstaller`

**Build:**
```bash
python build_app.py
```

**Or manually:**
```bash
pyinstaller --onefile --windowed --name=Primerly app.py
```

### B. Using the Desktop App Wrapper

**Run directly:**
```bash
python desktop_app.py
```

This creates a native desktop application with GUI that starts the Flask server and opens it in your browser.

## Option 3: Local Network Deployment

**Make accessible to other devices:**
```bash
# Edit app.py to change host
app.run(host='0.0.0.0', port=5000, debug=False)
```

## Production Considerations

### Security
- Change the secret key in `app.py`
- Use HTTPS in production
- Set up proper firewall rules
- Regular security updates

### Performance
- Use Gunicorn for production
- Set up caching (Redis)
- Optimize database queries
- Monitor resource usage

### Monitoring
- Set up logging
- Monitor application health
- Track user analytics
- Error reporting

## Troubleshooting

### Common Issues
1. **Port already in use**: Change port in `app.py`
2. **Dependencies missing**: Run `pip install -r requirements.txt`
3. **Permission errors**: Check file permissions
4. **Import errors**: Verify Python path

### Debug Mode
For development, enable debug mode:
```python
app.run(debug=True)
```

## Support

For deployment issues:
- Contact: jeevithacm21@gmail.com
- Website: https://www.biovagon.org/
- GitHub Issues: Create an issue on the repository

---

**Author**: Jeevitha C M, Founder of Biovagon  
**Organization**: Biovagon  
**Website**: https://www.biovagon.org/  
**License**: MIT 