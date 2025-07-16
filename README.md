# HAISim (Hospital-Acquired infection Interventions Simulator) Shiny App

This repository contains the code for **HAISim**, a Shiny app designed to simulate interventions for hospital-acquired infections.

## 🏥 Overview  

HAISim allows users to run simulations focusing on:      

- **Enhanced treatment intervention ONLY**  
- **Implementation of both Enhanced Treatment AND Infection Prevention Interventions**  

These simulations help assess the impacts on key outcomes such as:  
✔ Mortality  
✔ Lives saved  
✔ Changes in length of stay  
✔ Reductions in patient-days 


## 📥 Installation  

1. Download all the files in this repository.  
2. Open **`run.R`** in your R environment.  
3. Install the required packages (lines 31–37 in `run.R`).  
4. Run the entire script and specify the directory where the app is located.  

## 📂 Files  

- **`run.R`** – Main script to run the app.  
- **`app.R`** – Loads packages, sources other files, and initializes the UI.  
- **`functions.R`** – Contains all essential functions used in the app.  
- **`interface.R`** – Defines the app's user interface.  
- **`www/`** – Directory containing graphics and other static assets.  

## 🛠️ License  

**Copyright © 2025 Jean-Pierre Gnimatin, Marlon Grodd, Susanne Weber, Derek Hazard, Martin Wolkewitz**  

This program is free software: you can redistribute and/or modify it under the terms of the **GNU General Public License (GPL),** either version 3 or any later version.  

This program is distributed in the hope that it will be useful, but without any warranty, without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.  

You should have received a copy of the **GNU General Public License** along with this program. If not, visit [www.gnu.org/licenses](http://www.gnu.org/licenses/).  
