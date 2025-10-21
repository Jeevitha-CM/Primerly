#!/usr/bin/env python3
"""
Minimal Flask app for testing Render deployment
"""
from flask import Flask, render_template, jsonify
import os

app = Flask(__name__)
app.secret_key = 'test-secret-key'

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/health')
def health():
    return jsonify({'status': 'healthy', 'message': 'Primerly is working!'})

@app.route('/test')
def test():
    return jsonify({'message': 'Test endpoint working', 'version': '1.0.0'})

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    debug_mode = os.environ.get("FLASK_ENV") != "production"
    app.run(host='0.0.0.0', port=port, debug=debug_mode)
