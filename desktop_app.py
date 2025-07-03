import tkinter as tk
from tkinter import ttk, messagebox
import webbrowser
import threading
import time
from app import app

class PrimerlyDesktopApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Primerly - PCR Primer Design Tool")
        self.root.geometry("400x300")
        self.root.resizable(False, False)
        
        # Center the window
        self.center_window()
        
        # Create GUI
        self.create_widgets()
        
        # Flask app thread
        self.flask_thread = None
        self.server_running = False
        
    def center_window(self):
        """Center the window on screen"""
        self.root.update_idletasks()
        width = self.root.winfo_width()
        height = self.root.winfo_height()
        x = (self.root.winfo_screenwidth() // 2) - (width // 2)
        y = (self.root.winfo_screenheight() // 2) - (height // 2)
        self.root.geometry(f'{width}x{height}+{x}+{y}')
        
    def create_widgets(self):
        """Create the main GUI widgets"""
        # Main frame
        main_frame = ttk.Frame(self.root, padding="20")
        main_frame.grid(row=0, column=0, sticky="nsew")
        
        # Title
        title_label = ttk.Label(main_frame, text="Primerly", font=("Arial", 24, "bold"))
        title_label.grid(row=0, column=0, columnspan=2, pady=(0, 10))
        
        subtitle_label = ttk.Label(main_frame, text="PCR Primer Design Tool", font=("Arial", 12))
        subtitle_label.grid(row=1, column=0, columnspan=2, pady=(0, 30))
        
        # Status
        self.status_label = ttk.Label(main_frame, text="Ready to start", font=("Arial", 10))
        self.status_label.grid(row=2, column=0, columnspan=2, pady=(0, 20))
        
        # Buttons
        self.start_button = ttk.Button(main_frame, text="Start Primerly", command=self.start_server)
        self.start_button.grid(row=3, column=0, padx=(0, 10), pady=10)
        
        self.stop_button = ttk.Button(main_frame, text="Stop Server", command=self.stop_server, state="disabled")
        self.stop_button.grid(row=3, column=1, padx=(10, 0), pady=10)
        
        self.open_browser_button = ttk.Button(main_frame, text="Open in Browser", command=self.open_browser, state="disabled")
        self.open_browser_button.grid(row=4, column=0, columnspan=2, pady=10)
        
        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress.grid(row=5, column=0, columnspan=2, sticky="ew", pady=10)
        
    def start_server(self):
        """Start the Flask server in a separate thread"""
        if not self.server_running:
            self.flask_thread = threading.Thread(target=self.run_flask, daemon=True)
            self.flask_thread.start()
            
            self.start_button.config(state="disabled")
            self.stop_button.config(state="normal")
            self.status_label.config(text="Starting server...")
            self.progress.start()
            
            # Check if server is running
            self.root.after(1000, self.check_server_status)
    
    def run_flask(self):
        """Run Flask app"""
        try:
            app.run(host='127.0.0.1', port=5000, debug=False, use_reloader=False)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to start server: {str(e)}")
            self.stop_server()
    
    def check_server_status(self):
        """Check if the server is running"""
        try:
            import requests
            response = requests.get("http://127.0.0.1:5000", timeout=1)
            if response.status_code == 200:
                self.server_running = True
                self.status_label.config(text="Server running on http://127.0.0.1:5000")
                self.progress.stop()
                self.open_browser_button.config(state="normal")
            else:
                self.root.after(1000, self.check_server_status)
        except:
            self.root.after(1000, self.check_server_status)
    
    def stop_server(self):
        """Stop the Flask server"""
        self.server_running = False
        self.start_button.config(state="normal")
        self.stop_button.config(state="disabled")
        self.open_browser_button.config(state="disabled")
        self.status_label.config(text="Server stopped")
        self.progress.stop()
    
    def open_browser(self):
        """Open the app in default browser"""
        webbrowser.open("http://127.0.0.1:5000")

def main():
    root = tk.Tk()
    app = PrimerlyDesktopApp(root)
    root.mainloop()

if __name__ == "__main__":
    main() 