#!/usr/bin/env python
# coding: utf-8

# In[1]:


pip install matplotlib numpy


# In[ ]:


import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# --- Physics Logic & Calculations ---

def calculate_external_flow(density, viscosity, conductivity, prandtl, velocity, length, temp_fluid, temp_surface):
    """
    Calculates h for External Flow over a Flat Plate.
    Logic: Determines regime based on Reynolds number (Re_L).
    """
    # Calculate Reynolds Number at the end of the plate
    Re_L = (density * velocity * length) / viscosity
    
    # Determine Flow Regime and Nusselt Number (Nu)
    # Transition typically assumed at Re = 5e5
    if Re_L < 500000:
        regime = "Laminar"
        # Correlation for Laminar flow over flat plate
        Nu_L = 0.664 * (Re_L**0.5) * (prandtl**(1/3))
    else:
        regime = "Turbulent (Mixed)"
        # Correlation for Mixed Boundary Layer
        Nu_L = (0.037 * (Re_L**0.8) - 871) * (prandtl**(1/3))
    
    # Calculate Heat Transfer Coefficient (h)
    h = (Nu_L * conductivity) / length

    # --- Data for Graphing ---
    # Generates data for plotting Temperature vs. Length
    x_vals = np.linspace(0, length, 100)
    
    return h, Re_L, regime, x_vals

def calculate_internal_flow(density, viscosity, conductivity, prandtl, velocity, diameter, temp_inlet, temp_surface, length):
    """
    Calculates h for Internal Flow inside a circular pipe (Constant Surface Temp).
    Logic: Determines regime based on Reynolds number (Re_d).
    """
    # Calculate Reynolds Number
    Re_d = (density * velocity * diameter) / viscosity
    
    # Determine Flow Regime
    # Transition for internal flow is typically Re = 2300
    if Re_d < 2300:
        regime = "Laminar"
        Nu_d = 3.66 # Constant for fully developed laminar flow (Constant Ts)
    else:
        regime = "Turbulent"
        # Dittus-Boelter Equation (assuming heating for n=0.4)
        n = 0.4 if temp_surface > temp_inlet else 0.3
        Nu_d = 0.023 * (Re_d**0.8) * (prandtl**n)

    # Calculate Heat Transfer Coefficient (h)
    h = (Nu_d * conductivity) / diameter

    # --- Data for Graphing (Temperature Distribution) ---
    # We calculate the Mean Fluid Temperature (Tm) as it increases along the pipe
    # First, derive Specific Heat (Cp) from Prandtl number: Pr = Cp * mu / k
    cp = (prandtl * conductivity) / viscosity
    
    # Mass flow rate calculation
    area = np.pi * (diameter / 2)**2
    m_dot = density * velocity * area
    perimeter = np.pi * diameter
    
    x_vals = np.linspace(0, length, 100)
    # Equation for Tm(x) with constant surface temperature
    tm_vals = temp_surface - (temp_surface - temp_inlet) * np.exp(-(perimeter * h * x_vals) / (m_dot * cp))
    
    return h, Re_d, regime, x_vals, tm_vals

# --- GUI Application ---

class HeatTransferApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Convection Heat Transfer Calculator (MECH 3310)")
        self.root.geometry("1000x750")

        # 1. Input Frame
        input_frame = ttk.LabelFrame(root, text="Input Parameters", padding="15")
        input_frame.pack(side="left", fill="y", padx=15, pady=15)

        # Selection for External vs Internal
        ttk.Label(input_frame, text="Select Convection Type:", font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky="w", pady=10)
        self.mode_var = tk.StringVar(value="External (Flat Plate)")
        mode_cb = ttk.Combobox(input_frame, textvariable=self.mode_var, values=("External (Flat Plate)", "Internal (Pipe)"), state="readonly")
        mode_cb.grid(row=0, column=1, pady=10)
        mode_cb.bind("<<ComboboxSelected>>", self.update_inputs)

        # Dynamic Input Fields Area
        self.entries = {}
        self.param_frame = ttk.Frame(input_frame)
        self.param_frame.grid(row=1, column=0, columnspan=2, sticky="nsew")
        
        self.create_input_fields()

        # Calculate Button
        calc_btn = ttk.Button(input_frame, text="CALCULATE & PLOT", command=self.calculate)
        calc_btn.grid(row=2, column=0, columnspan=2, pady=25, sticky="ew")

        # Results Display Text
        self.result_text = tk.StringVar()
        result_label = ttk.Label(input_frame, textvariable=self.result_text, justify="left", font=("Consolas", 11), background="#f0f0f0", relief="sunken", padding=10)
        result_label.grid(row=3, column=0, columnspan=2, sticky="nsew")

        # 2. Graph Frame
        self.graph_frame = ttk.Frame(root)
        self.graph_frame.pack(side="right", fill="both", expand=True, padx=15, pady=15)
        
        self.fig, self.ax = plt.subplots(figsize=(5, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.graph_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def create_input_fields(self):
        """Generates input fields based on the selected mode."""
        for widget in self.param_frame.winfo_children():
            widget.destroy()
        self.entries.clear()

        # Standard Fluid Properties
        common_params = [
            ("Density (kg/m^3)", "1.225"),
            ("Dyn. Viscosity (kg/m.s)", "1.8e-5"),
            ("Thermal Cond. (W/m.K)", "0.025"),
            ("Prandtl Number", "0.7"),
            ("Velocity (m/s)", "5.0"),
        ]

        # Mode-Specific Parameters
        if self.mode_var.get() == "External (Flat Plate)":
            specific_params = [
                ("Plate Length (m)", "1.0"),
                ("Fluid Temp (C)", "25"),
                ("Surface Temp (C)", "80")
            ]
        else: # Internal
            specific_params = [
                ("Pipe Diameter (m)", "0.05"),
                ("Pipe Length (m)", "2.0"),
                ("Inlet Temp (C)", "25"),
                ("Surface Temp (C)", "80")
            ]

        # Draw the fields
        row = 0
        for label_text, default_val in common_params + specific_params:
            ttk.Label(self.param_frame, text=label_text).grid(row=row, column=0, sticky="w", pady=5)
            entry = ttk.Entry(self.param_frame)
            entry.insert(0, default_val)
            entry.grid(row=row, column=1, pady=5)
            # Store entry widget in dictionary for easy access
            key = label_text.split(" (")[0]
            self.entries[key] = entry
            row += 1

    def update_inputs(self, event):
        self.create_input_fields()

    def calculate(self):
        try:
            # 1. Retrieve Values
            rho = float(self.entries["Density"].get())
            mu = float(self.entries["Dyn. Viscosity"].get())
            k = float(self.entries["Thermal Cond."].get())
            pr = float(self.entries["Prandtl Number"].get())
            vel = float(self.entries["Velocity"].get())
            ts = float(self.entries["Surface Temp"].get())

            # 2. Perform Calculations based on Mode
            if "External" in self.mode_var.get():
                L = float(self.entries["Plate Length"].get())
                tf = float(self.entries["Fluid Temp"].get())
                
                h, Re, regime, x_vals = calculate_external_flow(rho, mu, k, pr, vel, L, tf, ts)
                
                # Format Results
                res_str = (f"--- RESULTS ---\n\n"
                           f"Flow Regime: {regime}\n"
                           f"Reynolds No: {Re:.2e}\n\n"
                           f"Heat Transfer Coeff (h):\n"
                           f"{h:.2f} W/m^2.K")
                self.result_text.set(res_str)

                # Plot External (Surface vs Fluid Temp)
                self.ax.clear()
                self.ax.plot(x_vals, [ts]*len(x_vals), label="Surface Temp (Ts)", color='red', linewidth=2)
                self.ax.plot(x_vals, [tf]*len(x_vals), label="Fluid Temp (T_inf)", color='blue', linestyle="--")
                self.ax.set_title("Temperature Profile (External Flow)")
                self.ax.set_xlabel("Distance along Plate (m)")
                self.ax.set_ylabel("Temperature (°C)")
                self.ax.legend()
                self.ax.grid(True, linestyle='--', alpha=0.7)

            else: # Internal
                D = float(self.entries["Pipe Diameter"].get())
                L = float(self.entries["Pipe Length"].get())
                tin = float(self.entries["Inlet Temp"].get())

                h, Re, regime, x_vals, tm_vals = calculate_internal_flow(rho, mu, k, pr, vel, D, tin, ts, L)

                # Format Results
                res_str = (f"--- RESULTS ---\n\n"
                           f"Flow Regime: {regime}\n"
                           f"Reynolds No: {Re:.2e}\n\n"
                           f"Heat Transfer Coeff (h):\n"
                           f"{h:.2f} W/m^2.K")
                self.result_text.set(res_str)

                # Plot Internal (Mean Fluid Temp Rising)
                self.ax.clear()
                self.ax.plot(x_vals, tm_vals, label="Mean Fluid Temp (Tm)", color='blue', linewidth=2)
                self.ax.plot(x_vals, [ts]*len(x_vals), label="Surface Temp (Ts)", color='red', linestyle="--")
                self.ax.set_title("Temperature Profile (Internal Flow)")
                self.ax.set_xlabel("Distance along Pipe (m)")
                self.ax.set_ylabel("Temperature (°C)")
                self.ax.legend()
                self.ax.grid(True, linestyle='--', alpha=0.7)

            self.canvas.draw()

        except ValueError:
            messagebox.showerror("Input Error", "Please check your inputs. Ensure all fields contain numbers.")

# --- Main Execution ---
if __name__ == "__main__":
    root = tk.Tk()
    app = HeatTransferApp(root)
    root.mainloop()


# In[ ]:




