import tkinter as tk
from tkinter import messagebox

from code_files.pilot_cancer import single_chromosome_process


def open_no_window():
    # Create a new window
    no_window = tk.Tk()
    no_window.title("Generate Interval Tables and Plots")

    # Create fields for user inputs
    tk.Label(no_window, text="Input File Path:").grid(row=0, column=0)
    input_file_entry = tk.Entry(no_window)
    input_file_entry.grid(row=0, column=1)

    tk.Label(no_window, text="Reference Type:").grid(row=1, column=0)
    reference_type_entry = tk.Entry(no_window)
    reference_type_entry.grid(row=1, column=1)

    tk.Label(no_window, text="Save Directory:").grid(row=2, column=0)
    save_directory_entry = tk.Entry(no_window)
    save_directory_entry.grid(row=2, column=1)

    tk.Label(no_window, text="Invert (Fill if True):").grid(row=3, column=0)
    invert_entry = tk.Entry(no_window)
    invert_entry.grid(row=3, column=1)

    tk.Label(no_window, text="Window Size:").grid(row=4, column=0)
    window_size_entry = tk.Entry(no_window)
    window_size_entry.grid(row=4, column=1)

    tk.Label(no_window, text="Error Size:").grid(row=5, column=0)
    error_size_entry = tk.Entry(no_window)
    error_size_entry.grid(row=5, column=1)

    # Function to handle button click
    def execute_function():
        input_file = input_file_entry.get()
        reference_type = reference_type_entry.get()
        save_directory = save_directory_entry.get()
        invert = invert_entry.get().lower() in ['true', '1', 't', 'y', 'yes']
        window_size = int(window_size_entry.get())
        error_size = int(error_size_entry.get())

        try:
            # Assume create_tables_and_plots is already imported and available
            from code_files.pilot_cancer import create_tables_and_plots
            create_tables_and_plots(input_file, reference_type, save_directory, invert, window_size, error_size)
            messagebox.showinfo("Success", "Operation completed successfully!")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # Button to initiate the process
    execute_button = tk.Button(no_window, text="Execute", command=execute_function)
    execute_button.grid(row=6, columnspan=2)

    no_window.mainloop()


def open_yes_window():
    # Create a new window
    yes_window = tk.Tk()
    yes_window.title("Generate Interval Tables and Plots")

    # Create fields for user inputs
    tk.Label(yes_window, text="Input File Path:").grid(row=0, column=0)
    input_file_entry = tk.Entry(yes_window)
    input_file_entry.grid(row=0, column=1)

    tk.Label(yes_window, text="Reference Type:").grid(row=1, column=0)
    reference_type_entry = tk.Entry(yes_window)
    reference_type_entry.grid(row=1, column=1)

    tk.Label(yes_window, text="Save Directory:").grid(row=2, column=0)
    save_directory_entry = tk.Entry(yes_window)
    save_directory_entry.grid(row=2, column=1)

    tk.Label(yes_window, text="Save Directory Plots:").grid(row=3, column=0)
    save_directory_plot_entry = tk.Entry(yes_window)
    save_directory_plot_entry.grid(row=3, column=1)

    tk.Label(yes_window, text="Invert (Fill if True):").grid(row=4, column=0)
    invert_entry = tk.Entry(yes_window)
    invert_entry.grid(row=4, column=1)

    tk.Label(yes_window, text="Chromosome Number:").grid(row=5, column=0)
    chrom_number_entry = tk.Entry(yes_window)
    chrom_number_entry.grid(row=5, column=1)

    tk.Label(yes_window, text="Window Size:").grid(row=6, column=0)
    window_size_entry = tk.Entry(yes_window)
    window_size_entry.grid(row=6, column=1)

    tk.Label(yes_window, text="Error Size:").grid(row=7, column=0)
    error_size_entry = tk.Entry(yes_window)
    error_size_entry.grid(row=7, column=1)

    # Function to handle button click
    def execute_function():
        input_file = input_file_entry.get()
        reference_type = reference_type_entry.get()
        save_directory = save_directory_entry.get()
        save_directory_plot = save_directory_plot_entry.get()
        invert = invert_entry.get().lower() in ['true', '1', 't', 'y', 'yes']
        chrom_number = chrom_number_entry.get()
        window_size = int(window_size_entry.get())
        error_size = int(error_size_entry.get())

        try:
            # Assume create_tables_and_plots is already imported and available
            from code_files.pilot_cancer import create_tables_and_plots
            single_chromosome_process(input_file, reference_type, save_directory, save_directory_plot, invert, chrom_number, window_size, error_size)
            messagebox.showinfo("Success", "Operation completed successfully!")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    # Button to initiate the process
    execute_button = tk.Button(yes_window, text="Execute", command=execute_function)
    execute_button.grid(row=8, columnspan=2)

    yes_window.mainloop()


def handle_response(is_yes):
    if is_yes:
        open_yes_window()
    else:
        open_no_window()


# Create the main window
root = tk.Tk()
root.title("Chromosome Selection")
root.geometry("500x300")

# Add a label with the question
label = tk.Label(root, text="Do you want to work on a specific chromosome?")
label.pack(pady=20)

# Add 'Yes' button
yes_button = tk.Button(root, text="Yes", command=lambda: handle_response(True))
yes_button.pack(side=tk.LEFT, padx=20, pady=20)

# Add 'No' button
no_button = tk.Button(root, text="No", command=lambda: handle_response(False))
no_button.pack(side=tk.RIGHT, padx=20, pady=20)

# Start the Tkinter event loop
root.mainloop()
