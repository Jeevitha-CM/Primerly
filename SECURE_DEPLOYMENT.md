# Secure Deployment Guide for Primerly

## 🔒 Files Protected from GitHub

The following sensitive files are **NEVER** committed to GitHub:

### User Data & Analytics
- `site_stats.json` - Page views, visitor counts, order interest
- `early_access.csv` - User signup data (names, emails, phones)
- `user_data/` - Any user-generated content
- `analytics/` - Usage analytics and tracking data

### Security & Configuration
- `*.key`, `*.pem`, `*.p12`, `*.pfx` - SSL certificates and private keys
- `secrets.json`, `config.json` - API keys and configuration
- `.env.production`, `.env.local` - Environment variables
- Any file containing passwords, API keys, or tokens

## 📁 Files Safe to Deploy

These files are **safe** to commit to GitHub:

### Core Application
- `app.py` - Main Flask application
- `wsgi.py` - WSGI entry point
- `requirements.txt` - Python dependencies
- `Procfile` - Heroku deployment config
- `render.yaml` - Render deployment config

### Templates & Static
- `templates/` - HTML templates
- `static/` - CSS, JS, images (except sensitive logos)
- `README.md` - Documentation
- `LICENSE` - License file

### Configuration
- `.gitignore` - Git ignore rules
- `DEPLOYMENT_GUIDE.md` - Deployment instructions
- `FEATURES.md` - Feature documentation

## 🚀 Deployment Steps

### 1. Clean Your Repository
```bash
# Remove sensitive files from Git tracking
git rm --cached site_stats.json early_access.csv
git rm --cached *.key *.pem secrets.json

# Add to .gitignore (already done)
echo "site_stats.json" >> .gitignore
echo "early_access.csv" >> .gitignore
```

### 2. Environment Variables for Production
Set these in your hosting platform (Render/Heroku):

```bash
# Required
FLASK_ENV=production
SECRET_KEY=your-super-secret-key-here

# Optional
PORT=5000
NCBI_EMAIL=your-email@domain.com
```

### 3. Deploy to Render
1. Connect your GitHub repo
2. Set environment variables in Render dashboard
3. Deploy automatically on push

### 4. Deploy to Heroku
```bash
# Create app
heroku create your-app-name

# Set environment variables
heroku config:set FLASK_ENV=production
heroku config:set SECRET_KEY=your-secret-key

# Deploy
git push heroku main
```

## 🔐 Security Best Practices

### 1. Never Commit Sensitive Data
- Use environment variables for all secrets
- Store user data in databases, not files
- Use `.env` files locally (gitignored)

### 2. Production Configuration
```python
# In app.py - use environment variables
import os
app.secret_key = os.environ.get('SECRET_KEY', 'dev-key-only')
Entrez.email = os.environ.get('NCBI_EMAIL', 'default@example.com')
```

### 3. Database for Production
Replace file-based storage with proper database:
- PostgreSQL (Render/Heroku)
- SQLite (local development)
- Store user data, stats, early access signups

### 4. HTTPS Only
- Force HTTPS in production
- Use secure cookies
- Validate all inputs

## 📊 Data Handling

### Local Development
- `site_stats.json` and `early_access.csv` are created locally
- These files are gitignored and won't be committed

### Production
- Use database tables instead of JSON/CSV files
- Implement proper user data protection
- Follow GDPR/privacy regulations

## 🚨 What NOT to Deploy

❌ **Never commit:**
- User personal information
- Analytics data
- API keys or secrets
- Database files
- Log files with sensitive data
- SSL certificates
- Configuration with passwords

✅ **Always deploy:**
- Source code
- Templates and static assets
- Documentation
- Configuration files (without secrets)
- Dependencies list

## 🔄 After Deployment

1. **Verify sensitive files are not in repository:**
   ```bash
   git ls-files | grep -E "\.(json|csv|key|pem)$"
   ```

2. **Test production environment:**
   - Check all features work
   - Verify no sensitive data is exposed
   - Test user data collection

3. **Monitor and maintain:**
   - Regular security updates
   - Monitor for data leaks
   - Backup user data securely

---

**Remember:** When in doubt, don't commit it. It's better to be overly cautious with sensitive data.
