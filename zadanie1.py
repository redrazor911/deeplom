import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import os
import sys
from skbio.io import read
import pandas as pd
from Bio import Entrez
import argparse
import time
from datetime import datetime, timedelta

def get_virus_taxids(names_filename, nodes_filename, parent_taxons=None, parent_taxons_file=None):
    with open(names_filename, 'r', encoding='utf-8') as names_file:
        try:
            names_df = read(names_file, format='taxdump', into=pd.DataFrame, scheme="names")
        except Exception as e:
            raise ValueError(f"Error reading names.dmp: {e}")



    with open(nodes_filename, 'r', encoding='utf-8') as nodes_file:
        try:
            nodes_df = read(nodes_file, format='taxdump', into=pd.DataFrame, scheme="nodes_slim")
        except Exception as e:
            raise ValueError(f"Error reading nodes.dmp: {e}")

  
    names_df = names_df[names_df['name_class'] == 'scientific name']


    if parent_taxons_file:
        try:
            with open(parent_taxons_file, 'r', encoding='utf-8') as f:
                file_taxons = [line.strip() for line in f]
                if parent_taxons:
                    parent_taxons.extend(file_taxons)
                else:
                    parent_taxons = file_taxons

        except FileNotFoundError:
            raise FileNotFoundError(f"Parent taxons file not found: {parent_taxons_file}")

    if not parent_taxons:
        parent_taxons = ['Viruses']  

    virus_taxids = []
    for taxid in nodes_df.index:
        try:
            lineage = get_lineage(taxid, nodes_df)
            lineage_names = names_df.loc[lineage, 'name_txt'].tolist()

          
            if any(taxon in lineage_names for taxon in parent_taxons):
                virus_taxids.append(taxid)
        except KeyError:
            print(f"TaxID {taxid} not found in names_df.  Skipping.")


    return sorted(virus_taxids)

def get_lineage(taxid, nodes_df):
    lineage = []
    current_id = taxid
    while current_id != 1: 
        if current_id not in nodes_df.index:
            break
        lineage.append(current_id)
        current_id = nodes_df.loc[current_id, 'parent_tax_id']
    lineage.append(1) 
    return lineage


def browse_file(entry):
   
    filename = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, filename)


def run_analysis():
  
    names_filename = names_entry.get()
    nodes_filename = nodes_entry.get()
    output_taxids_filename = output_taxids_entry.get()
    parent_taxons_file = parent_taxons_file_entry.get()
    parent_taxons_text = parent_taxons_text_area.get("1.0", tk.END).strip()   
    parent_taxons = [taxon.strip() for taxon in parent_taxons_text.split('\n') if taxon.strip()]



    status_text.delete("1.0", tk.END) 
    try:
        status_text.insert(tk.END, "Starting analysis...\n")
        window.update()

    
        start_time = time.time()
        status_text.insert(tk.END, "Estimating analysis time...\n") 
        window.update()

     
        sample_size = min(100, len(pd.read_csv(nodes_filename, sep='\t', header=None, low_memory=False)))  
        sample_nodes_df = pd.read_csv(nodes_filename, sep='\t', header=None, nrows=sample_size, low_memory=False)
        sample_taxids = sample_nodes_df.iloc[:, 0].tolist() 
        num_taxids = len(sample_taxids)
        estimated_time_per_taxid = 0
        start_sample_time = time.time()
        taxids = get_virus_taxids(names_filename, nodes_filename, parent_taxons, parent_taxons_file)
        end_sample_time = time.time()
        estimated_time_per_taxid = (end_sample_time - start_sample_time) / num_taxids 

        total_nodes_df = pd.read_csv(nodes_filename, sep='\t', header=None, low_memory=False)
        num_all_taxids = len(total_nodes_df)

        estimated_total_time = estimated_time_per_taxid * num_all_taxids
        estimated_time_str = str(timedelta(seconds=estimated_total_time))  
        status_text.insert(tk.END, f"Estimated analysis time: {estimated_time_str}\n")
        window.update()

        taxids = get_virus_taxids(names_filename, nodes_filename, parent_taxons, parent_taxons_file)

        if output_taxids_filename:
            with open(output_taxids_filename, "w") as f:
                for taxid in taxids:
                    f.write(str(taxid) + "\n")
            status_text.insert(tk.END, f"Virus TaxIDs written to {output_taxids_filename}\n")
            window.update()



        status_text.insert(tk.END, "Analysis complete!\n") 
        window.update()

        end_time = time.time()
        actual_time = end_time - start_time 
        actual_time_str = str(timedelta(seconds=actual_time)) 
        status_text.insert(tk.END, f"Actual analysis time: {actual_time_str}\n") 
        window.update()

    except FileNotFoundError as e:
        messagebox.showerror("Error", str(e))
        status_text.insert(tk.END, str(e) + "\n")
    except ValueError as e:
        messagebox.showerror("Error", str(e))
        status_text.insert(tk.END, str(e) + "\n")
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred: {e}")
        status_text.insert(tk.END, f"An unexpected error occurred: {e}\n")
    finally:
        window.update()


window = tk.Tk()
window.title("Virus TaxID Extractor") 


names_label = tk.Label(window, text="names.dmp:")
names_label.grid(row=0, column=0, sticky=tk.W)
names_entry = tk.Entry(window, width=50)
names_entry.grid(row=0, column=1, sticky=(tk.W + tk.E))
names_browse_button = tk.Button(window, text="Browse", command=lambda: browse_file(names_entry))
names_browse_button.grid(row=0, column=2, sticky=tk.W)

nodes_label = tk.Label(window, text="nodes.dmp:")
nodes_label.grid(row=1, column=0, sticky=tk.W)
nodes_entry = tk.Entry(window, width=50)
nodes_entry.grid(row=1, column=1, sticky=(tk.W + tk.E))
nodes_browse_button = tk.Button(window, text="Browse", command=lambda: browse_file(nodes_entry))
nodes_browse_button.grid(row=1, column=2, sticky=tk.W)

output_taxids_label = tk.Label(window, text="Output TaxIDs File:")
output_taxids_label.grid(row=2, column=0, sticky=tk.W)
output_taxids_entry = tk.Entry(window, width=50)
output_taxids_entry.grid(row=2, column=1, sticky=(tk.W + tk.E))
output_taxids_browse_button = tk.Button(window, text="Browse", command=lambda: browse_file(output_taxids_entry))
output_taxids_browse_button.grid(row=2, column=2, sticky=tk.W)



parent_taxons_file_label = tk.Label(window, text="Parent Taxons File:")
parent_taxons_file_label.grid(row=3, column=0, sticky=tk.W)
parent_taxons_file_entry = tk.Entry(window, width=50)
parent_taxons_file_entry.grid(row=3, column=1, sticky=(tk.W + tk.E))
parent_taxons_file_browse_button = tk.Button(window, text="Browse", command=lambda: browse_file(parent_taxons_file_entry))
parent_taxons_file_browse_button.grid(row=3, column=2, sticky=tk.W)

parent_taxons_text_label = tk.Label(window, text="Parent Taxons (one per line):")
parent_taxons_text_label.grid(row=4, column=0, sticky=tk.W)
parent_taxons_text_area = scrolledtext.ScrolledText(window, width=50, height=5)
parent_taxons_text_area.grid(row=4, column=1, sticky=(tk.W+tk.E))


run_button = tk.Button(window, text="Run Analysis", command=run_analysis)
run_button.grid(row=5, column=1, pady=10)


status_label = tk.Label(window, text="Status:")
status_label.grid(row=6, column=0, sticky=tk.W)
status_text = scrolledtext.ScrolledText(window, width=60, height=10)
status_text.grid(row=6, column=1, columnspan=2, sticky=(tk.W + tk.E + tk.N + tk.S))



window.columnconfigure(1, weight=1)
window.rowconfigure(6, weight=1)


window.mainloop()
