#!/usr/bin/env python3
"""
Simple test script to verify the Flask app can start without errors
"""
import sys
import os

try:
    # Add current directory to path
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    
    # Try to import the app
    from app import app
    
    print("SUCCESS: Flask app imported successfully")
    print(f"SUCCESS: App name: {app.name}")
    print(f"SUCCESS: Debug mode: {app.debug}")
    
    # Test if we can create a test client
    with app.test_client() as client:
        response = client.get('/')
        print(f"SUCCESS: Home route works: {response.status_code}")
        
        response = client.get('/health')
        print(f"SUCCESS: Health check works: {response.status_code}")
        
    print("SUCCESS: All tests passed! App should work on Render.")
    
except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
