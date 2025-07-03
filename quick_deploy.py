#!/usr/bin/env python3
"""
Quick Deployment Script for Primerly
This script helps you deploy Primerly to free hosting services for public access.
"""

import os
import sys
import subprocess
import webbrowser
from pathlib import Path

def print_header():
    print("🚀 Primerly Quick Deployment")
    print("=" * 50)
    print("Make your primer design tool accessible to everyone!")
    print()

def check_git():
    """Check if git is initialized and ready"""
    try:
        result = subprocess.run(['git', 'status'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ Git repository is ready")
            return True
        else:
            print("❌ Git repository not initialized")
            return False
    except FileNotFoundError:
        print("❌ Git not found. Please install Git first.")
        return False

def check_files():
    """Check if all required files exist"""
    required_files = ['app.py', 'requirements.txt', 'wsgi.py', 'Procfile']
    missing_files = []
    
    for file in required_files:
        if not Path(file).exists():
            missing_files.append(file)
    
    if missing_files:
        print(f"❌ Missing required files: {', '.join(missing_files)}")
        return False
    else:
        print("✅ All required files present")
        return True

def show_deployment_options():
    """Show deployment options"""
    print("\n📋 Choose your deployment option:")
    print()
    print("1. 🎯 Render.com (Recommended - Free & Easy)")
    print("   - No credit card required")
    print("   - Automatic deployments")
    print("   - Custom domain support")
    print()
    print("2. 🚄 Railway.app (Fast & Simple)")
    print("   - Very fast deployment")
    print("   - Auto-detects Flask")
    print("   - Free tier available")
    print()
    print("3. 🐘 Heroku (Popular Choice)")
    print("   - Well-established platform")
    print("   - Good documentation")
    print("   - Free tier available")
    print()
    print("4. 🐍 PythonAnywhere (Python-focused)")
    print("   - Python-specific hosting")
    print("   - Easy setup")
    print("   - Free tier available")
    print()

def open_deployment_guide():
    """Open deployment guide in browser"""
    guide_path = Path("DEPLOYMENT_SUMMARY.md")
    if guide_path.exists():
        print("\n📖 Opening detailed deployment guide...")
        try:
            webbrowser.open(f"file://{guide_path.absolute()}")
        except:
            print("Could not open guide automatically. Please open DEPLOYMENT_SUMMARY.md manually.")
    else:
        print("❌ Deployment guide not found")

def prepare_for_deployment():
    """Prepare files for deployment"""
    print("\n🔧 Preparing for deployment...")
    
    # Check if .gitignore exists
    if not Path('.gitignore').exists():
        print("Creating .gitignore...")
        gitignore_content = """# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
env/
venv/
.venv/
ENV/
env.bak/
venv.bak/

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db

# Logs
*.log

# Environment variables
.env
"""
        with open('.gitignore', 'w') as f:
            f.write(gitignore_content)
    
    # Check if Procfile exists
    if not Path('Procfile').exists():
        print("Creating Procfile...")
        with open('Procfile', 'w') as f:
            f.write('web: gunicorn app:app')
    
    print("✅ Deployment files ready!")

def main():
    print_header()
    
    # Check prerequisites
    if not check_git():
        print("\n💡 To initialize git repository:")
        print("   git init")
        print("   git add .")
        print("   git commit -m 'Initial commit'")
        return
    
    if not check_files():
        print("\n❌ Please ensure all required files are present before deploying.")
        return
    
    prepare_for_deployment()
    show_deployment_options()
    
    print("\n🎯 RECOMMENDED: Start with Render.com")
    print("   It's the easiest option for beginners!")
    print()
    
    choice = input("Enter your choice (1-4) or 'q' to quit: ").strip()
    
    if choice == '1':
        print("\n🚀 Opening Render.com...")
        webbrowser.open("https://render.com")
        print("1. Sign up for a free account")
        print("2. Click 'New Web Service'")
        print("3. Connect your GitHub repository")
        print("4. Configure:")
        print("   - Build Command: pip install -r requirements.txt")
        print("   - Start Command: gunicorn app:app")
        print("5. Deploy and get your public URL!")
        
    elif choice == '2':
        print("\n🚄 Opening Railway.app...")
        webbrowser.open("https://railway.app")
        print("1. Sign up for a free account")
        print("2. Click 'New Project'")
        print("3. Connect your GitHub repository")
        print("4. Railway will auto-detect Flask and deploy!")
        
    elif choice == '3':
        print("\n🐘 Opening Heroku...")
        webbrowser.open("https://heroku.com")
        print("1. Sign up for a free account")
        print("2. Install Heroku CLI")
        print("3. Run: heroku create primerly-app")
        print("4. Run: git push heroku main")
        
    elif choice == '4':
        print("\n🐍 Opening PythonAnywhere...")
        webbrowser.open("https://pythonanywhere.com")
        print("1. Sign up for a free account")
        print("2. Upload your Primerly files")
        print("3. Configure WSGI file to point to app.py")
        print("4. Get your public URL!")
        
    elif choice.lower() == 'q':
        print("Goodbye! 👋")
        return
        
    else:
        print("Invalid choice. Please try again.")
        return
    
    print("\n📚 Need more help? Check the detailed deployment guide!")
    open_deployment_guide()
    
    print("\n🎉 Once deployed, you'll get a public URL that you can share with friends!")
    print("   Example: https://primerly.onrender.com")

if __name__ == "__main__":
    main() 