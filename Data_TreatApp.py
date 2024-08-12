import pandas as pd
import matplotlib.pyplot as plt
import tkinter as tk
import scipy.signal
from scipy.signal import find_peaks, welch
from scipy.integrate import solve_ivp, trapezoid, cumulative_trapezoid, quad, quad_vec
from scipy.interpolate import CubicSpline
import scipy.fftpack
import numpy as np
import rainflow
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from tkinter import Tk, messagebox, filedialog, StringVar, IntVar, Checkbutton, Button, OptionMenu, Canvas, Frame, Scrollbar, ttk
from tkinter import simpledialog as sd
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from pandas.plotting import register_matplotlib_converters
from apread import APReader
import os
import folium
import sys
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QDialog
from PyQt5.QtCore import QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView

register_matplotlib_converters()

# class HTMLViewer(QDialog):
#     def __init__(self, html_file_path):
#         super().__init__()

#         web_view = QWebEngineView()
#         web_view.setUrl(QUrl.fromLocalFile(os.path.abspath(html_file_path)))

#         # Create a layout and set it as central
#         layout = QVBoxLayout()
#         layout.addWidget(web_view)

#         central_widget = QWidget()
#         central_widget.setLayout(layout)

        # self.setLayout(layout)
        
class DataVisualizationApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Data Visualization App")
        
        #HTMLviewer tracker
        self.html_viewer = None  # Add this line to keep track of the HTML viewer window
        
        # DataFrame
        self.df = pd.DataFrame()
        self.file_name = 'None'
        
        #some cubicsplines to be
        self.a_rd_func = None
        self.a_rt_func = None
        
        #Parametros gerais do sistema moto, peso piloto, peso carga, entre eixos
        self.sys_param = None
        
        # Initialize canvas attribute
        self.canvas = None  # Add this line
        
        # Create main frame
        main_frame = Frame(root)
        main_frame.pack(fill='both', expand=True)

        # Create left canvas for buttons
        left_canvas = Canvas(main_frame, width=root.winfo_screenwidth() * 0.3, height=root.winfo_screenheight() * 0.9)
        left_canvas.grid(row=0, column=0, sticky="nsw")

        # Load Data Button
        self.load_data_button = Button(left_canvas, text="Load Data", command=self.load_data)
        self.load_data_button.grid(row=0, column=0, padx=10, pady=10, sticky="nw")

        # X-Axis Dropdown Menu
        self.x_axis_var = StringVar(root)
        self.x_axis_var.set("Select X-Axis")
        self.x_axis_dropdown = OptionMenu(left_canvas, self.x_axis_var, "Select X-Axis")
        self.x_axis_dropdown.grid(row=1, column=0, padx=10, pady=10, sticky="w")

        # Checkboxes
        self.checkbox_vars = []
        self.checkboxes = []  # List to store checkbox widgets

        # Scrollbar and Canvas for Checkboxes
        scrollbar_height_percentage = 0.6 # Adjust this value as needed
        scrollbar_height = int(scrollbar_height_percentage * root.winfo_screenheight())

        self.canvas_frame = Frame(left_canvas)
        self.scrollbar = Scrollbar(self.canvas_frame, orient="vertical")
        self.checkbox_canvas = Canvas(self.canvas_frame, yscrollcommand=self.scrollbar.set, height=scrollbar_height)
        self.checkbox_canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")
        self.canvas_frame.grid(row=2, column=0, padx=10, pady=10, sticky="w")

        # Plot Selected Data Button
        self.plot_button = Button(left_canvas, text="Plot Selected Data", command=self.plot_data)
        self.plot_button.grid(row=3, column=0, padx=10, pady=10, sticky="nw")
      
        # Botão Otimização de parâmetros de suspensão
        self.plot_vibes_button = Button(left_canvas, text = "Otimizar param suspensão", command = self.start_suspParamOptim)
        self.plot_vibes_button.grid(row=3, column=0, padx=10, pady=10, sticky="ne")
        
        # Botão Vibes process data 
        self.plot_vibes_button = Button(left_canvas, text = "Resposta 2GDL", command = self.vibes_process_data)
        self.plot_vibes_button.grid(row=4, column=0, padx=10, pady=10, sticky="ne")
        
        # Slice Data Button
        self.slice_data_button = Button(left_canvas, text="Slice Data", command=self.pre_slice_data)
        self.slice_data_button.grid(row=6, column=0, padx=10, pady=10, sticky="ne")
        
        #Rainflow Count Button
        #self.filter_data_Button = Button(left_canvas, text="Rainflow Count", command=self.Rainflow_count)
        #self.filter_data_Button.grid(row=4, column=0, padx=10, pady=10, sticky="nw")
        
        #PSD Button
        self.filter_data_Button = Button(left_canvas, text="PSD", command=self.PSD)
        self.filter_data_Button.grid(row=4, column=0, padx=10, pady=10, sticky="nw")
        
        #ISO 2631 Calc Button
        self.filter_data_Button = Button(left_canvas, text="Param. ISO 2631", command=self.ISO_2631)
        self.filter_data_Button.grid(row=5, column=0, padx=10, pady=10, sticky="nw")
        
        #Histogram Plot Button
        #self.filter_data_Button = Button(left_canvas, text="Histogram Plot", command=self.Histogram_plot)
        #self.filter_data_Button.grid(row=6, column=0, padx=10, pady=10, sticky="nw")
        
        # Create checkboxes after canvas_frame is defined
        self.create_checkboxes()

        # Create right canvas for plots
        right_canvas = Canvas(main_frame, width=root.winfo_screenwidth() * 0.6, height=root.winfo_screenheight() * 0.9)
        right_canvas.grid(row=0, column=1, sticky="nsew")

        # Plot Frame
        plot_frame_height_percentage = 0.9  # Adjust this value as needed
        plot_frame_width_percentage = 0.6  # Adjust this value as needed
        plot_frame_height = int(plot_frame_height_percentage * root.winfo_screenheight())
        plot_frame_width = int(plot_frame_width_percentage * root.winfo_screenwidth())

        self.plot_frame = Frame(right_canvas, height=plot_frame_height, width=plot_frame_width)
        self.plot_frame.pack(padx=10, pady=10, fill='both', expand=True)
        self.index_zoomed_subplot = 0
        self.axes_subplot = None
        self.previous_ax = None
        self.click_count = 0
        self.syncro_zooms = 0
        self.t_ini = 0
        self.t_end = 0
        
    def load_data(self):
        file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("BIN files", "*.bin")])
        print(file_path)
        if file_path:
            
            file_ext = os.path.splitext(file_path)[1].lower()
            self.file_name = os.path.splitext(os.path.basename(file_path))[0]
            
            if file_ext == '.csv':
                self.df = pd.read_csv(file_path)
                #self.process_data()
                self.df.columns = self.df.columns.str.strip()
                self.populate_dropdown()
                self.create_checkboxes()
            elif file_ext == '.bin':
            
                file = APReader(file_path)
                
                for chx in file.Channels:
                    
                    self.df[chx.Name] = chx.data
                self.df.columns = self.df.columns.str.strip()
                #self.process_data()
                self.populate_dropdown()
                self.create_checkboxes()
                # Perform any additional processing for .BIN files if needed
            else:
                print("Unsupported file type. Please select a .CSV or .BIN file.")
         
    def process_data(self):
        # 3.1 and 3.2: Calculate derivatives
        self.df['Vel_Diant'] = self.df['Desl_Diant'].diff() / self.df['Time  2 - default sample rate'].diff()
        self.df['Vel_Tras'] = self.df['Desl_Tras'].diff() / self.df['Time  2 - default sample rate'].diff()
        self.df['Ac_Ang_x'] = self.df['Rotation_A_Rotation_X'].diff()/self.df['Time  2 - default sample rate'].diff()
        self.df['Ac_Ang_y'] = self.df['Rotation_A_Rotation_Y'].diff()/self.df['Time  2 - default sample rate'].diff()
        #3.3: Convert GPS minutes to decimals of degrees
        self.df['GPS_Latitude'] = -1 * (self.df['GPS_PositionLatitude_GPS_Latitude_Degree'] + (
                self.df['GPS_PositionLatitude_GPS_Latitude_Minutes'] / 60.0))
        self.df['GPS_Longitude'] = -1 * (self.df['GPS_PositionLongitude_GPS_Longitude_Degree'] + (
                self.df['GPS_PositionLongitude_GPS_Longitude_Minutes'] / 60.0))
        
        #self.df['Force_Diant'] = self.df['Desl_Diant']*self.stiff_k_bengala + self.const_a_bengala*self.df['Vel_Diant']**(self.const_b_bengala)
        #self.df['Force_Tras'] = self.df['Desl_Tras']*self.stiff_k_cofap + self.const_a_cofap*self.df['Vel_Trans']**(self.const_b_cofap)
        
        Total_Run_Distance = ((self.df['GPS_CourseSpeed_GPS_Speed']/3.6)*self.df['Time  2 - default sample rate'].diff()).sum()
        self.df['Displacement'] = ((self.df['GPS_CourseSpeed_GPS_Speed']/3.6)*self.df['Time  2 - default sample rate'].diff()).cumsum()
        
        print('Distância de coleta = ', Total_Run_Distance,'m')
        
        
        #Coeficientes para converter deformação em carga para o chassi v2.5
        m =  1.829332830976119  
        b =  5.749999999999971
        
        #Convertendo deformação e em carga
        self.df['Carga eixo dianteiro'] = (self.df['G14b'])*m + b
        #self.slice_data()
    
    def find_Ixx_cg_coord_forOptim(self,m_piloto,m_carga,entre_eixos):
        # Define the arrays for Ppil and Pcarga
        Ppil_values = np.array([65, 75, 85, 95, 105, 115])
        Pcarga_values = np.array([5, 10, 15, 20, 25])

        # Define your 6x5 matrices of calculated values for the Ixx, b, and h parameters
        Ixx_matrix = np.array([
            [26.035286371, 27.93070260, 29.70560258, 31.37161510, 32.93891942],
            [27.53764919, 29.37375382, 31.10054961, 32.72794644, 34.26469143],
            [28.99200350, 30.77741687, 32.46290809, 34.05701972, 35.567347575],
            [30.40679880, 32.14842653, 33.79809577, 35.36324433, 36.85052900],
            [31.78861167, 33.49206497, 35.11039254, 36.65012728, 38.11715031],
            [33.14263741, 34.81253335, 36.40322425, 37.92049279, 39.36957211]
        ])

        b_matrix = np.array([
            [636.87, 620.45, 604.85, 590.01, 575.88],
            [628.00, 612.59, 597.93, 583.94, 570.59],
            [619.97, 605.47, 591.64, 578.41, 565.77],
            [612.68, 598.99, 585.89, 573.36, 561.35],
            [606.02, 593.05, 580.63, 568.71, 557.28],
            [599.92, 587.61, 575.79, 564.44, 553.52]
        ])

        h_matrix = np.array([
            [691.91, 700.57, 708.79, 716.62, 724.07],
            [705.89, 713.78, 721.30, 728.47, 735.31],
            [718.54, 725.77, 732.67, 739.27, 745.58],
            [730.03, 736.69, 743.05, 749.14, 754.98],
            [740.53, 746.67, 752.56, 758.21, 763.64],
            [750.15, 755.85, 761.32, 766.57, 771.62]
        ])

        # Function to find the two closest values in an array to a given value
        def find_closest_values(array, value):
            idx = np.searchsorted(array, value)
            if idx == 0:
                return array[0], array[1]
            if idx == len(array):
                return array[-2], array[-1]
            return array[idx - 1], array[idx]

        # Function to calculate the proportional rule for a value between two points
        def proportional_rule(value, x1, x2, y1, y2):
            return y1 + (value - x1) * (y2 - y1) / (x2 - x1)

        # Function to calculate the parameter value based on input values and parameter matrix
        def calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, parameter_matrix):
            # Find the two closest values of Ppil and Pcarga in the arrays
            closest_Ppil_1, closest_Ppil_2 = find_closest_values(Ppil_values, Ppil_input)
            closest_Pcarga_1, closest_Pcarga_2 = find_closest_values(Pcarga_values, Pcarga_input)

            # Find the corresponding parameter values in the parameter_matrix
            param_11 = parameter_matrix[np.where(Ppil_values == closest_Ppil_1), np.where(Pcarga_values == closest_Pcarga_1)][0][0]
            param_12 = parameter_matrix[np.where(Ppil_values == closest_Ppil_1), np.where(Pcarga_values == closest_Pcarga_2)][0][0]
            param_21 = parameter_matrix[np.where(Ppil_values == closest_Ppil_2), np.where(Pcarga_values == closest_Pcarga_1)][0][0]
            param_22 = parameter_matrix[np.where(Ppil_values == closest_Ppil_2), np.where(Pcarga_values == closest_Pcarga_2)][0][0]

            # Calculate the proportional rule for the given Ppil and Pcarga
            param_value = proportional_rule(Ppil_input, closest_Ppil_1, closest_Ppil_2,
                                            proportional_rule(Pcarga_input, closest_Pcarga_1, closest_Pcarga_2, param_11, param_12),
                                            proportional_rule(Pcarga_input, closest_Pcarga_1, closest_Pcarga_2, param_21, param_22))
            
            return param_value

        # Prompt the user to input Ppil and Pcarga values
        # Ppil_input = float(input("Enter the value of Ppil: "))
        # Pcarga_input = float(input("Enter the value of Pcarga: "))
        # Wb_input = float(input("Enter the value of Wheelbase: "))
        
                    
        
        Ppil_input = m_piloto
        Pcarga_input = m_carga
        Wb_input = entre_eixos
        
        # Calculate the values for Ixx, b, and h
        Ixx_value = calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, Ixx_matrix)
        b_value = calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, b_matrix)
        h_value = calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, h_matrix)
        mass_value = Ppil_input+Pcarga_input
        p_value = Wb_input-b_value
        mass_rear = mass_value*(b_value/Wb_input)
        mass_front = mass_value*((Wb_input-b_value)/Wb_input)
            
 
        self.sys_param = [Ixx_value, b_value, h_value, mass_value, p_value,mass_rear,mass_front]
            
    def find_Ixx_cg_coord(self):
        # Define the arrays for Ppil and Pcarga
        Ppil_values = np.array([65, 75, 85, 95, 105, 115])
        Pcarga_values = np.array([5, 10, 15, 20, 25])

        # Define your 6x5 matrices of calculated values for the Ixx, b, and h parameters
        Ixx_matrix = np.array([
            [26.035286371, 27.93070260, 29.70560258, 31.37161510, 32.93891942],
            [27.53764919, 29.37375382, 31.10054961, 32.72794644, 34.26469143],
            [28.99200350, 30.77741687, 32.46290809, 34.05701972, 35.567347575],
            [30.40679880, 32.14842653, 33.79809577, 35.36324433, 36.85052900],
            [31.78861167, 33.49206497, 35.11039254, 36.65012728, 38.11715031],
            [33.14263741, 34.81253335, 36.40322425, 37.92049279, 39.36957211]
        ])

        b_matrix = np.array([
            [636.87, 620.45, 604.85, 590.01, 575.88],
            [628.00, 612.59, 597.93, 583.94, 570.59],
            [619.97, 605.47, 591.64, 578.41, 565.77],
            [612.68, 598.99, 585.89, 573.36, 561.35],
            [606.02, 593.05, 580.63, 568.71, 557.28],
            [599.92, 587.61, 575.79, 564.44, 553.52]
        ])

        h_matrix = np.array([
            [691.91, 700.57, 708.79, 716.62, 724.07],
            [705.89, 713.78, 721.30, 728.47, 735.31],
            [718.54, 725.77, 732.67, 739.27, 745.58],
            [730.03, 736.69, 743.05, 749.14, 754.98],
            [740.53, 746.67, 752.56, 758.21, 763.64],
            [750.15, 755.85, 761.32, 766.57, 771.62]
        ])

        # Function to find the two closest values in an array to a given value
        def find_closest_values(array, value):
            idx = np.searchsorted(array, value)
            if idx == 0:
                return array[0], array[1]
            if idx == len(array):
                return array[-2], array[-1]
            return array[idx - 1], array[idx]

        # Function to calculate the proportional rule for a value between two points
        def proportional_rule(value, x1, x2, y1, y2):
            return y1 + (value - x1) * (y2 - y1) / (x2 - x1)

        # Function to calculate the parameter value based on input values and parameter matrix
        def calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, parameter_matrix):
            # Find the two closest values of Ppil and Pcarga in the arrays
            closest_Ppil_1, closest_Ppil_2 = find_closest_values(Ppil_values, Ppil_input)
            closest_Pcarga_1, closest_Pcarga_2 = find_closest_values(Pcarga_values, Pcarga_input)

            # Find the corresponding parameter values in the parameter_matrix
            param_11 = parameter_matrix[np.where(Ppil_values == closest_Ppil_1), np.where(Pcarga_values == closest_Pcarga_1)][0][0]
            param_12 = parameter_matrix[np.where(Ppil_values == closest_Ppil_1), np.where(Pcarga_values == closest_Pcarga_2)][0][0]
            param_21 = parameter_matrix[np.where(Ppil_values == closest_Ppil_2), np.where(Pcarga_values == closest_Pcarga_1)][0][0]
            param_22 = parameter_matrix[np.where(Ppil_values == closest_Ppil_2), np.where(Pcarga_values == closest_Pcarga_2)][0][0]

            # Calculate the proportional rule for the given Ppil and Pcarga
            param_value = proportional_rule(Ppil_input, closest_Ppil_1, closest_Ppil_2,
                                            proportional_rule(Pcarga_input, closest_Pcarga_1, closest_Pcarga_2, param_11, param_12),
                                            proportional_rule(Pcarga_input, closest_Pcarga_1, closest_Pcarga_2, param_21, param_22))
            
            return param_value

        # Prompt the user to input Ppil and Pcarga values
        # Ppil_input = float(input("Enter the value of Ppil: "))
        # Pcarga_input = float(input("Enter the value of Pcarga: "))
        # Wb_input = float(input("Enter the value of Wheelbase: "))
        
        # New window for user inputs
        control_var = tk.IntVar() #Variável de controle que faz a rotina aguardar o input do usuário
        input_window = tk.Toplevel(root)
        input_window.title("Input de massas e entre-eixos")
        
        # Create a frame for the input fields and button in the new window
        frame = ttk.Frame(input_window, padding="10")
        frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Field 1
        ttk.Label(frame, text="Massa do piloto [kg]:").grid(row=0, column=0, sticky=tk.W)
        entry1 = ttk.Entry(frame, width=20)
        entry1.grid(row=0, column=1)

        # Field 2
        ttk.Label(frame, text="Massa da carga [kg]:").grid(row=1, column=0, sticky=tk.W)
        entry2 = ttk.Entry(frame, width=20)
        entry2.grid(row=1, column=1)
        
        # Field 3
        ttk.Label(frame, text="Entre eixos [mm]:").grid(row=2, column=0, sticky=tk.W)
        entry3 = ttk.Entry(frame, width=20)
        entry3.grid(row=2, column=1)
                    
        def handle_submit():
            Ppil_input = float(entry1.get())
            Pcarga_input = float(entry2.get())
            Wb_input = float(entry3.get())
            # Calculate the values for Ixx, b, and h
            Ixx_value = calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, Ixx_matrix)
            b_value = calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, b_matrix)
            h_value = calculate_parameter(Ppil_input, Pcarga_input, Ppil_values, Pcarga_values, h_matrix)
            mass_value = Ppil_input + Pcarga_input + 125
            p_value = Wb_input-b_value
            mass_rear = mass_value*(b_value/Wb_input)
            mass_front = mass_value*((Wb_input-b_value)/Wb_input)
            
            print("The calculated value of Ixx:", Ixx_value)
            print("The calculated value of b:", b_value)
            print("The calculated value of p:", p_value)
            print("The calculated value of h:", h_value)
            print("The calculated value of mass:", mass_value)
            print('Massa roda dianteira:', mass_front)
            print('Massa na roda traseira:', mass_rear)
            self.sys_param = [Ixx_value, b_value, h_value, mass_value, p_value,mass_rear,mass_front]
            
            control_var.set(1)
            input_window.destroy()
                
        # Submit Button
        submit_button = ttk.Button(frame, text="Submit", command=handle_submit)
        submit_button.grid(row=3, column=0, columnspan=2)
        
        # Wait for the user to press submit
        input_window.wait_variable(control_var)
        
    def vibes_process_data(self):
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return    
        
        self.find_Ixx_cg_coord()
        #Parâmetros invariantes do sistema:
        m_total = self.sys_param[3] #[kg]
        I_total = self.sys_param[0] # [kg*m^2]
        #messagebox.showinfo("Resultado", f"Momento de inércia total I_tot [kg*m^2] = {I_total: .2f}")
        
        b_value = self.sys_param[1]/1000 # distância centro de massa até eixo traseiro[mm]
        h_value = self.sys_param[2]/1000 # altura centro de massa com relação ao solo [mm]
        p_value = self.sys_param[4]/1000 # distância centro de massa até o eixo dianteiro [mm]
        m_f = self.sys_param[6] # massa apoiada na roda da frente [kg]
        m_r = self.sys_param[5] # massa apoiada na roda de trás [kg]
        
        # Parâmetros geométricos da suspensão, inserir a baixo 
        #Dianteira
        caster = np.radians(28)
        h_diant = 0.5
        x_diant = h_diant*np.sin(caster)
        
        #Traseira
        phi = np.radians(20.18)
        veta_ini = np.radians(172.11)
        L1 = 301.12/1000
        L = 393/1000
        h_do = 48.66 /1000
        Lsa_ini = 325.7/1000
        r_do = np.sqrt(h_do**2 + L1**2)
        gama_do = np.arctan(h_do/L1)
        x_ma = -r_do*np.cos(veta_ini + gama_do) - Lsa_ini*np.sin(phi)
        y_ma = Lsa_ini*np.cos(phi) - r_do*np.sin(veta_ini + gama_do)
        
        
        def pondera_acc(acc,freq):
            #Coeficientes para um polinomio de 10° Graus ajustado para representar tabela 3 de Freq vs Pesos da norma ISO2631
            Wk = pd.DataFrame()
            # Valores retirados da ISO2631 Tabela  3
            Wk['Frequências'] = [ 0.1, 0.125, 0.16, 0.2, 0.25, 0.315, 0.4, 0.5, 0.63, 0.8,   1, 1.25, 1.6,   2, 2.5, 3.15,   4,    5,  6.3,    8,  10, 12.5,  16,  20,  25, 31.5,  40,  50,  63,  80,  100, 125,  160,  200, 250,  315,  400]
            Wk['Pesos'] =       [  10,  10,   10, 10,  10,   10, 10, 10,  10, 10, 10,  10, 10, 10, 631,  804, 967, 1039, 1054, 1036, 988,  902, 768, 636, 513,  405, 314, 246, 186, 132, 88.7,  54, 28.5, 15.2, 7.9, 3.98, 1.95]
            for i in range(len(Wk['Pesos'])):
                if Wk.loc[i,'Frequências'] < 2:
                    Wk.loc[i,'Pesos'] = 0
                elif 2 <= Wk.loc[i,'Frequências'] <= 8:
                    Wk.loc[i, 'Pesos'] = 1000
                else:
                    Wk.loc[i, 'Pesos'] = 1000
            
            Wk_freq = Wk['Frequências'].values
            Wk_pesos = (Wk['Pesos'].values)/1000
            cs = CubicSpline(Wk_freq, Wk_pesos)
            peso = cs(freq)

            acc_w = acc*peso

            return acc_w        

        def pondera_acc_ISO2631(acc,freq):
            #Coeficientes para um polinomio de 10° Graus ajustado para representar tabela 3 de Freq vs Pesos da norma ISO2631
            Wk = pd.DataFrame()
            # Valores retirados da ISO2631 Tabela  3
            Wk['Frequências'] = [ 0.1, 0.125, 0.16, 0.2, 0.25, 0.315, 0.4, 0.5, 0.63, 0.8,   1, 1.25, 1.6,   2, 2.5, 3.15,   4,    5,  6.3,    8,  10, 12.5,  16,  20,  25, 31.5,  40,  50,  63,  80,  100, 125,  160,  200, 250,  315,  400]
            Wk['Pesos'] =       [31.2,  48.6,   79, 121,  182,   263, 352, 418,  459, 477, 482,  484, 494, 531, 631,  804, 967, 1039, 1054, 1036, 988,  902, 768, 636, 513,  405, 314, 246, 186, 132, 88.7,  54, 28.5, 15.2, 7.9, 3.98, 1.95]
            Wk_freq = Wk['Frequências'].values
            Wk_pesos = (Wk['Pesos'].values)/1000
            cs = CubicSpline(Wk_freq, Wk_pesos)
            peso = cs(freq)

            acc_w = acc*peso

            return acc_w   
        
        def freq_weighting(y,t):

            fft_y_complete = scipy.fft.fft(y)
        
            N = fft_y_complete.size
        
            freq_y_complete = scipy.fft.fftfreq(N,t)
        
            fft_y = fft_y_complete[1:N//2]
            freq_y = freq_y_complete[1:N//2]
            fft_y_w = np.copy(fft_y_complete)
            fft_y_w[1:N//2] = pondera_acc(fft_y, freq_y)

            if N%2 == 0:
                                
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2][::-1])
            else:
                
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2+1][::-1])
                
            y_w = np.real(scipy.fft.ifft(fft_y_w))        
            return y_w        
        
        def freq_weighting_ISO2631(y,t):

            fft_y_complete = scipy.fft.fft(y)
        
            N = fft_y_complete.size
        
            freq_y_complete = scipy.fft.fftfreq(N,t)
        
            fft_y = fft_y_complete[1:N//2]
            freq_y = freq_y_complete[1:N//2]
            fft_y_w = np.copy(fft_y_complete)
            fft_y_w[1:N//2] = pondera_acc_ISO2631(fft_y, freq_y)

            if N%2 == 0:
                                
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2][::-1])
            else:
                
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2+1][::-1])
                
            y_w = np.real(scipy.fft.ifft(fft_y_w))        
            return y_w 
        
        #Definindo excitações ao longo do tempo
        # Get the selected X-axis column from the dropdown menu
        x_axis_column = self.x_axis_var.get()
        t_span = (self.df[x_axis_column].min(), self.df[x_axis_column].max()) #Intervalo de tempo da coleta
        t_eval = self.df[x_axis_column].values #Array contendo as leituras de tempo
        t = 1/600 #Período de aquisição 1/f[Hz]
        
        acc_rodadiant = sd.askstring('Input','Digite o nome da série de acelerações da roda dianteira: ')
        acc_rodatras = sd.askstring('Input','Digite o nome da série de acelerações da roda traseira: ')
        

        
        acc_rd = (self.df[acc_rodadiant]*9.81).values #Array contendo as leituras de aceleração na roda dianteira
        acc_rd = acc_rd - acc_rd.mean()
              
        acc_rt = (self.df[acc_rodatras]*9.81).values #Array contendo as leituras de aceleração na roda traseira
        acc_rt = acc_rt - acc_rt.mean()
        
        if 'Ac_Amortecedor-z' in self.df.columns:
            acc_ma = (self.df['Ac_Amortecedor-z']*9.81).values
            accma_w = freq_weighting_ISO2631(acc_ma,t)
        else:
            print('Não foram carregados dados de acelerações no montante do amortecedor')
            
        # Pré definindo vetores para receber os resultados das integrais das acelerações
        vel_rd = np.zeros_like(acc_rd)
        Desl_rd = np.zeros_like(acc_rd)
        vel_rt = np.zeros_like(acc_rt)
        Desl_rt = np.zeros_like(acc_rt)
        
        #Filtrando acelerações em baixas frequências
        accrd_w = freq_weighting(acc_rd, t)
        accrt_w = freq_weighting(acc_rt,t)
        
        
        vel_rd[1:] = cumulative_trapezoid(accrd_w, t_eval)
        vel_rt[1:] = cumulative_trapezoid(accrt_w,t_eval)
        
        vel_rd = vel_rd - vel_rd.mean()
        vel_rt = vel_rt - vel_rt.mean()
        vel_rd_cs = CubicSpline(t_eval, vel_rd)
        vel_rt_cs = CubicSpline(t_eval, vel_rt)
        
        Desl_rd[1:] = cumulative_trapezoid(vel_rd, t_eval)
        Desl_rt[1:] = cumulative_trapezoid(vel_rt, t_eval)
        Desl_rd_cs = CubicSpline(t_eval, Desl_rd)
        Desl_rt_cs = CubicSpline(t_eval, Desl_rt)
        
        # plt.plot(t_eval, Desl_rt, color = 'red')
        # plt.xlabel('time [s]')
        # plt.ylabel('deslocamento roda traseira [mm]')
        # plt.show()
        

         
        #Parâmetros variantes do sistema
        def rigidez_susp_frontal(desl_rd, z_rd, caster): #definindo rigidez equivalente  da suspensão frontal como função do deslocamento e dos parâmetros geométricos da suspensão
            #rigidez do componente "bengala"
            # desl_rel_f = desl_rd - z_rd
            
            # if desl_rel_f > 125:
            #     desl_rel_f = 125
            # elif desl_rel_f < 0:
            #     desl_rel_f = 0
                        
            # desl_susp_f = (desl_rel_f)/np.sin(caster)
            # coefs = [28928.73848161,  4779.11062105] # Coeficientes do polinomio que define a rigidez em função do deslocamento da suspensão
            # polynomial = np.poly1d(coefs)
            k_saf = 6400*2 #np.poly1d(polynomial)(desl_susp_f)
            k_f = k_saf/(np.cos(caster)**2)
         
            return k_f #valores dummy só pra rodar o código
        
        def rigidez_susp_traseira(desl_rt, z_rt, L1, L, Lsa, x_ma, y_ma, veta_ini): #definindo rigidez equivalente da suspensão traseira como função do deslocamento e dos parâmetros geométricos da suspensão
            
            #print('deslocamento roda traseira: ', desl_rt)
            veta_ = veta_ini + np.arcsin((desl_rt)/L)
            #print('angulo calculado pelo arcsen_rigidez: ', np.rad2deg(veta_))
            
            

            y_lsa = h_do*np.sin(veta_ - gama_do)
            x_lsa = np.sqrt(r_do**2 - y_lsa**2)
            
            desl_susp_rt = Lsa - np.sqrt((x_ma-x_lsa)**2 + (y_ma-y_lsa)**2)
           
            v_ratio = ((L1*(x_ma*np.sin(veta_)-y_ma*np.cos(veta_)))/(np.sqrt((L1**2)+(2*L1*((x_ma*np.cos(veta_))+(y_ma*np.sin(veta_))))+x_ma**2+y_ma**2)))*(1/(-L*np.cos(veta_)))
            
            coefs =  [ 1425.06336763, 21579.69218816, 62.2342865 ]
            polynomial = np.poly1d(coefs)
            k_sar = 21600*2 #np.polyder(polynomial)(desl_susp_rt)
            k_r = k_sar*v_ratio**2
            return k_r 
        
        def coef_amort_traseira(vel_rt, vz_rt, desl_rt, z_rt, L1, L, Lsa, x_ma, y_ma, veta_ini): #definindo amortecimento equivalente da suspensão traseira como função do deslocamento e dos parâmetros geométricos da suspensão
            
            
            veta_ = veta_ini + np.arcsin((desl_rt)/L)
            #print('angulo calculado pelo arcsen_amortecimento: ', np.rad2deg(veta_))
            y_lsa = h_do*np.sin(veta_ - gama_do)
            x_lsa = np.sqrt(r_do**2 - y_lsa**2)
            
            Lsa_i = np.sqrt((x_ma-x_lsa)**2 + (y_ma-y_lsa)**2)
            alpha = np.arccos((y_ma-y_lsa)/Lsa_i)
            v_susp_rt = ((vel_rt/np.cos(veta_))*L1/L)/np.cos(alpha)
           
            v_ratio = ((L1*(x_ma*np.sin(veta_)-y_ma*np.cos(veta_)))/(np.sqrt((L1**2)+(2*L1*((x_ma*np.cos(veta_))+(y_ma*np.sin(veta_))))+x_ma**2+y_ma**2)))*(1/(-L*np.cos(veta_)))
            coefs =  [  ]
            polynomial = np.poly1d(coefs)
            c_sar = 1350*2
            c_rt = c_sar*v_ratio**2
            return c_rt #valores dummy só pra rodar o código
                
        def coef_amort_frontal(vel_rd, vz_rd, caster): #definindo amortecimento equivalente da suspensão traseira como função do deslocamento e dos parâmetros geométricos da suspensão
            
            vel_susp_f = (vel_rd - vz_rd)/np.cos(caster)
            coefs =  [  ]
            polynomial = np.poly1d(coefs)
            c_saf = 2000*2 # Substituir pelo polinomio definidor do coeficiente de amortecidmento
            c_rd = c_saf/(np.cos(caster)**2)
            
            return c_rd #valores dummy só pra rodar o código
    
        #Definindo a funçao que carrega a equação do movimento
        def acc_2gdl(t, y):
            z, zdot, theta, thetadot = y
            
            z_rd = z - theta*p_value
            vz_rd = zdot - thetadot*p_value
            z_rt = z + theta*b_value
            vz_rt = zdot + thetadot*b_value
            # print('z_rd = ', z_rd, ' | vz_rd = ', vz_rd, ' | z_rt = ', z_rt, ' | vz_rt = ', vz_rt)
            
            desl_rd = Desl_rd_cs(t)
            desl_rt = Desl_rt_cs(t)
            Vel_rd = vel_rd_cs(t)
            Vel_rt = vel_rt_cs(t)
            # print('desl_rd = ', desl_rd,' | Vel_rd = ', Vel_rd, ' | desl_rt = ', desl_rt,' | Vel_rt = ', Vel_rt)
            
            kd = rigidez_susp_frontal(desl_rd,z_rd,caster)
            kt = rigidez_susp_traseira(desl_rt,z_rt,L1, L, Lsa_ini, x_ma, y_ma, veta_ini)
            cd = coef_amort_frontal(Vel_rd, vz_rd, caster)
            ct = coef_amort_traseira(Vel_rt, vz_rt, desl_rt, z_rt, L1, L, Lsa_ini, x_ma, y_ma, veta_ini)
            # print('kd = ', kd,' | kt = ', kt, ' | cd = ', cd,' | ct = ', ct)
            
            z_ddot = (-z*(kd+kt) - theta*(kt*b_value-kd*p_value) - zdot*(cd+ct) - thetadot*(ct*b_value-cd*p_value) + kd*desl_rd + kt*desl_rt + cd*Vel_rd + ct*Vel_rt - m_total*9.81) /m_total
            theta_ddot = (-z*(b_value*kt-p_value*kd) - theta*(kt*b_value**2+kd*p_value**2) - zdot*(b_value*ct-p_value*cd) - thetadot*(ct*b_value**2+cd*p_value**2) + kt*b_value*desl_rt - kd*p_value*desl_rd + ct*b_value*Vel_rt - cd*p_value*Vel_rd )/I_total
            # print('z_ddot = ', z_ddot, 'theta_ddot = ', theta_ddot)
               
            return z_ddot, theta_ddot
        
        y = [0,0,0,0]
        dt = 1/600
        z_array = np.zeros_like(t_eval)
        zdot_array = np.zeros_like(t_eval)
        theta_array = np.zeros_like(t_eval)
        thetadot_array = np.zeros_like(t_eval)
        zddot_array = np.zeros_like(t_eval)
        thetaddot_array = np.zeros_like(t_eval)
        acc_ma_array = np.zeros_like(t_eval)
        # count = 0
        # sol_2gdl = solve_ivp(acc_2gdl, t_span, y, args = (caster, L1 ,L , Lsa_ini, x_ma, y_ma, veta_ini, m_total, I_total, b_value, p_value, gama_do, h_do, r_do), t_eval = t_eval)
        # sol_z_dot = sol_2gdl.y[0]
        # sol_z_ddot = np.gradient(sol_z_dot,t_eval)
        # sol_theta_dot = sol_2gdl.y[1]
        # sol_theta_ddot = np.gradient(sol_theta_dot, t_eval)
        # sol_acc_ma = sol_z_ddot + sol_theta_ddot*b_value
        for i in tqdm(range(len(t_eval)), desc='Resolvendo sistema 2GDL'):
            z, zdot, theta, thetadot = y
            t = t_eval[i]
            
            z_ddot, theta_ddot = acc_2gdl(t,y)
            acc_ma = z_ddot + theta_ddot*b_value
            # if count >= 60000:
            #     messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            #     count = 0

            zdot += z_ddot*dt
            z += zdot*dt
            thetadot += theta_ddot*dt
            theta += thetadot*dt
            
            z_array[i] = z
            zdot_array[i] = zdot
            theta_array[i] = theta
            thetadot_array[i] = thetadot
            zddot_array[i] = z_ddot
            thetaddot_array[i] = theta_ddot
            acc_ma_array[i] = acc_ma
            y = z, zdot, theta, thetadot

        #     # count = count + 1
        # Create a figure and subplots
        plt.close('all')
        width_factor = 0.7
        fig = plt.figure(figsize=((root.winfo_screenwidth()/96)*width_factor, 6))

        plt.plot(t_eval, acc_ma_array, color = 'red',label = 'acc montante teorica')
        # plt.plot(t_eval, accma_w, color = 'blue',label = 'acc montante medida ponderada ISO')
        #plt.plot(t_eval, self.df['Ac_Amortecedor-z']*9.81, color = 'green',label = 'acc montante medida')
        plt.xlabel('time [s]')
        plt.ylabel('aceleração [m/s^2]')
        plt.legend()
        # Clear existing canvas and toolbar
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        if hasattr(self, "toolbar"):
            self.toolbar.destroy()

        # Create the canvas
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
            
        # Add navigation toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)    
        self.df['acc_z_ma_tr'] = acc_ma_array
               
    def start_suspParamOptim(self):
        
        self.find_Ixx_cg_coord()
        acc_rodadiant = sd.askstring('Input','Digite o nome da série de acelerações da roda dianteira: ')
        acc_rodatras = sd.askstring('Input','Digite o nome da série de acelerações da roda traseira: ')
        
        global progress, prog_bar
        progress = tk.Toplevel(root)
        progress.title('Iterando parametros de suspensão') 
        
        prog_bar = ttk.Progressbar(progress, orient='horizontal', length=300, mode='determinate')
        prog_bar.pack(pady=20)
        
        root.after(100, self.susp_paramOptim(acc_rodadiant, acc_rodatras))
               
    def susp_paramOptim(self,acc_r_diant,acc_r_tras):
        phi = np.radians(np.linspace(0.1,30,5)) #
        
        L1 = np.linspace(280,305,5)
        VDV = np.zeros((len(phi),len(L1)))
        smaller_VDV = 1000
        
        prog_bar['maximum'] = len(L1)*len(phi)
        prog_bar['value'] = 0
        

        
        for i in range(len(phi)):
            for j in range(len(L1)):
                print('Inclinação = ', np.degrees(phi[i]), 'L1 = ', L1[j])
                self.vibes_process_data_forOptim(phi[i], L1[j],acc_r_diant,acc_r_tras)
                VDV[i,j] = self.ISO_2631_forOptim()
                if VDV[i,j] < smaller_VDV:
                    smaller_VDV = VDV[i,j]
                    smaller_i = i
                    smaller_j = j
                print('\n VDV = ', VDV[i,j])
                
                prog_bar.step(1)
                progress.update()
                
        print('O menor valor de VDV é ', smaller_VDV, 'referente aos parâmetros, phi = ', np.degrees(phi[smaller_i]), 'L1 = ', L1[smaller_j])
        progress.destroy()
        
        plt.figure(figsize=(8,6))
        plt.imshow(VDV.T, extent=[np.degrees(phi.min()),np.degrees(phi.max()),L1.min(),L1.max()], origin = 'lower', aspect='auto', cmap='viridis')
        plt.colorbar(label='VDV')
        plt.xlabel('Ângulo do amortecedor[°]')
        plt.ylabel('Distância de ancoragem [mm]')
        plt.legend()
        plt.show()
                    
    def vibes_process_data_forOptim(self,phi,L_1,acc_rodadiant,acc_rodatras):
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return   
         

        
        
        #Parâmetros invariantes do sistema:
        m_total = self.sys_param[3] #[kg]
        I_total = self.sys_param[0] # [kg*m^2]
        #messagebox.showinfo("Resultado", f"Momento de inércia total I_tot [kg*m^2] = {I_total: .2f}")
        
        b_value = self.sys_param[1]/1000 # distância centro de massa até eixo traseiro[mm]
        h_value = self.sys_param[2]/1000 # altura centro de massa com relação ao solo [mm]
        p_value = self.sys_param[4]/1000 # distância centro de massa até o eixo dianteiro [mm]
        m_f = self.sys_param[6] # massa apoiada na roda da frente [kg]
        m_r = self.sys_param[5] # massa apoiada na roda de trás [kg]
        
        # Parâmetros geométricos da suspensão, inserir a baixo 
        #Dianteira
        caster = np.radians(28)
        h_diant = 0.5
        x_diant = h_diant*np.sin(caster)
        
        #Traseira
        veta_ini = np.radians(172.11)
        L1 = L_1/1000
        L = 393/1000
        h_do = 44.5 /1000
        Lsa_ini = 325.7/1000
        r_do = np.sqrt(h_do**2 + L1**2)
        gama_do = np.arctan(h_do/L1)
        x_ma = -r_do*np.cos(veta_ini + gama_do) - Lsa_ini*np.sin(phi)
        y_ma = Lsa_ini*np.cos(phi) - r_do*np.sin(veta_ini + gama_do)
        
        
        def pondera_acc(acc,freq):
            #Coeficientes para um polinomio de 10° Graus ajustado para representar tabela 3 de Freq vs Pesos da norma ISO2631
            Wk = pd.DataFrame()
            # Valores retirados da ISO2631 Tabela  3
            Wk['Frequências'] = [ 0.1, 0.125, 0.16, 0.2, 0.25, 0.315, 0.4, 0.5, 0.63, 0.8,   1, 1.25, 1.6,   2, 2.5, 3.15,   4,    5,  6.3,    8,  10, 12.5,  16,  20,  25, 31.5,  40,  50,  63,  80,  100, 125,  160,  200, 250,  315,  400]
            Wk['Pesos'] =       [  10,  10,   10, 10,  10,   10, 10, 10,  10, 10, 10,  10, 10, 10, 631,  804, 967, 1039, 1054, 1036, 988,  902, 768, 636, 513,  405, 314, 246, 186, 132, 88.7,  54, 28.5, 15.2, 7.9, 3.98, 1.95]
            for i in range(len(Wk['Pesos'])):
                if Wk.loc[i,'Frequências'] < 2:
                    Wk.loc[i,'Pesos'] = 0
                elif 2 <= Wk.loc[i,'Frequências'] <= 8:
                    Wk.loc[i, 'Pesos'] = 1000
                else:
                    Wk.loc[i, 'Pesos'] = 1000
            
            Wk_freq = Wk['Frequências'].values
            Wk_pesos = (Wk['Pesos'].values)/1000
            cs = CubicSpline(Wk_freq, Wk_pesos)
            peso = cs(freq)

            acc_w = acc*peso

            return acc_w        
        
        def freq_weighting(y,t):

            fft_y_complete = scipy.fft.fft(y)
        
            N = fft_y_complete.size
        
            freq_y_complete = scipy.fft.fftfreq(N,t)
        
            fft_y = fft_y_complete[1:N//2]
            freq_y = freq_y_complete[1:N//2]
            fft_y_w = np.copy(fft_y_complete)
            fft_y_w[1:N//2] = pondera_acc(fft_y, freq_y)

            if N%2 == 0:
                                
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2][::-1])
            else:
                
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2+1][::-1])
                
            y_w = np.real(scipy.fft.ifft(fft_y_w))        
            return y_w        
        
        #Definindo excitações ao longo do tempo
        # Get the selected X-axis column from the dropdown menu
        x_axis_column = self.x_axis_var.get()
        t_span = (self.df[x_axis_column].min(), self.df[x_axis_column].max()) #Intervalo de tempo da coleta
        t_eval = self.df[x_axis_column].values #Array contendo as leituras de tempo
        t = 1/600 #Período de aquisição 1/f[Hz]
        

        

        
        acc_rd = (self.df[acc_rodadiant]*9.81).values #Array contendo as leituras de aceleração na roda dianteira
        acc_rd = acc_rd - acc_rd.mean()
              
        acc_rt = (self.df[acc_rodatras]*9.81).values #Array contendo as leituras de aceleração na roda traseira
        acc_rt = acc_rt - acc_rt.mean()
        
        if 'Ac_Amortecedor-z' in self.df.columns:
            acc_ma = (self.df['Ac_Amortecedor-z']*9.81).values
            
        else:
            print('Aviso: Não foram carregados dados de acelerações no montante do amortecedor')
        
        # Pré definindo vetores para receber os resultados das integrais das acelerações
        vel_rd = np.zeros_like(acc_rd)
        Desl_rd = np.zeros_like(acc_rd)
        vel_rt = np.zeros_like(acc_rt)
        Desl_rt = np.zeros_like(acc_rt)
        
        #Filtrando acelerações em baixas frequências
        accrd_w = freq_weighting(acc_rd, t)
        accrt_w = freq_weighting(acc_rt,t)
        
        vel_rd[1:] = cumulative_trapezoid(accrd_w, t_eval)
        vel_rt[1:] = cumulative_trapezoid(accrt_w,t_eval)
        
        vel_rd = vel_rd - vel_rd.mean()
        vel_rt = vel_rt - vel_rt.mean()
        vel_rd_cs = CubicSpline(t_eval, vel_rd)
        vel_rt_cs = CubicSpline(t_eval, vel_rt)
        
        Desl_rd[1:] = cumulative_trapezoid(vel_rd, t_eval)
        Desl_rt[1:] = cumulative_trapezoid(vel_rt, t_eval)
        Desl_rd_cs = CubicSpline(t_eval, Desl_rd)
        Desl_rt_cs = CubicSpline(t_eval, Desl_rt)
       
        
        #Parâmetros variantes do sistema
        def rigidez_susp_frontal(desl_rd, z_rd, caster): #definindo rigidez equivalente  da suspensão frontal como função do deslocamento e dos parâmetros geométricos da suspensão
            #rigidez do componente "bengala"
            desl_rel_f = desl_rd - z_rd
            
            if desl_rel_f > 125:
                desl_rel_f = 125
            elif desl_rel_f < 0:
                desl_rel_f = 0
                        
            desl_susp_f = (desl_rel_f)/np.sin(caster)
            coefs = [28928.73848161,  4779.11062105] # Coeficientes do polinomio que define a rigidez em função do deslocamento da suspensão
            polynomial = np.poly1d(coefs)
            k_saf = 6400*2 #np.poly1d(polynomial)(desl_susp_f)
            k_f = k_saf/(np.cos(caster)**2)
         
            return k_f #valores dummy só pra rodar o código
        
        def rigidez_susp_traseira(desl_rt, z_rt, L1, L, Lsa, x_ma, y_ma, veta_ini): #definindo rigidez equivalente da suspensão traseira como função do deslocamento e dos parâmetros geométricos da suspensão
            
            #print('deslocamento roda traseira: ', desl_rt)
            veta_ = veta_ini + np.arcsin((desl_rt)/L)
            #print('angulo calculado pelo arcsen_rigidez: ', np.rad2deg(veta_))
            
            

            y_lsa = h_do*np.sin(veta_ - gama_do)
            x_lsa = np.sqrt(r_do**2 - y_lsa**2)
            
            desl_susp_rt = Lsa - np.sqrt((x_ma-x_lsa)**2 + (y_ma-y_lsa)**2)
           
            v_ratio = ((L1*(x_ma*np.sin(veta_)-y_ma*np.cos(veta_))))/(np.sqrt((L1**2)+(2*L1*((x_ma*np.cos((veta_))+(y_ma*np.sin(veta_))))+x_ma**2+y_ma**2)))*(1/(-L*np.cos((veta_))))
            
            coefs =  [ 1425.06336763, 21579.69218816, 62.2342865 ]
            polynomial = np.poly1d(coefs)
            k_sar = 21600*2 #np.polyder(polynomial)(desl_susp_rt)
            k_r = k_sar*v_ratio**2
            return k_r 
        
        def coef_amort_traseira(vel_rt, vz_rt, desl_rt, z_rt, L1, L, Lsa, x_ma, y_ma, veta_ini): #definindo amortecimento equivalente da suspensão traseira como função do deslocamento e dos parâmetros geométricos da suspensão
            
            
            veta_ = veta_ini + np.arcsin((desl_rt)/L)
            #print('angulo calculado pelo arcsen_amortecimento: ', np.rad2deg(veta_))
            y_lsa = h_do*np.sin(veta_ - gama_do)
            x_lsa = np.sqrt(r_do**2 - y_lsa**2)
            
            Lsa_i = np.sqrt((x_ma-x_lsa)**2 + (y_ma-y_lsa)**2)
            alpha = np.arccos((y_ma-y_lsa)/Lsa_i)
            v_susp_rt = ((vel_rt/np.cos(veta_))*L1/L)/np.cos(alpha)
           
            v_ratio = ((L1*(x_ma*np.sin(veta_)-y_ma*np.cos(veta_)))/(np.sqrt((L1**2)+(2*L1*((x_ma*np.cos(veta_))+(y_ma*np.sin(veta_))))+x_ma**2+y_ma**2)))*(1/(-L*np.cos(veta_)))
            coefs =  [  ]
            polynomial = np.poly1d(coefs)
            c_sar = 1350*2
            c_rt = c_sar*v_ratio**2
            return c_rt #valores dummy só pra rodar o código
                
        def coef_amort_frontal(vel_rd, vz_rd, caster): #definindo amortecimento equivalente da suspensão traseira como função do deslocamento e dos parâmetros geométricos da suspensão
            
            vel_susp_f = (vel_rd - vz_rd)/np.cos(caster)
            coefs =  [  ]
            polynomial = np.poly1d(coefs)
            c_saf = 2000*2 # Substituir pelo polinomio definidor do coeficiente de amortecidmento
            c_rd = c_saf/(np.cos(caster)**2)
            
            return c_rd #valores dummy só pra rodar o código
    
        #Definindo a funçao que carrega a equação do movimento
        def acc_2gdl(t, y):
            z, zdot, theta, thetadot = y
            
            z_rd = z - theta*p_value
            vz_rd = zdot - thetadot*p_value
            z_rt = z + theta*b_value
            vz_rt = zdot + thetadot*b_value
            # print('z_rd = ', z_rd, ' | vz_rd = ', vz_rd, ' | z_rt = ', z_rt, ' | vz_rt = ', vz_rt)
            
            desl_rd = Desl_rd_cs(t)
            desl_rt = Desl_rt_cs(t)
            Vel_rd = vel_rd_cs(t)
            Vel_rt = vel_rt_cs(t)
            # print('desl_rd = ', desl_rd,' | Vel_rd = ', Vel_rd, ' | desl_rt = ', desl_rt,' | Vel_rt = ', Vel_rt)
            
            kd = rigidez_susp_frontal(desl_rd,z_rd,caster)
            kt = rigidez_susp_traseira(desl_rt,z_rt,L1, L, Lsa_ini, x_ma, y_ma, veta_ini)
            cd = coef_amort_frontal(Vel_rd, vz_rd, caster)
            ct = coef_amort_traseira(Vel_rt, vz_rt, desl_rt, z_rt, L1, L, Lsa_ini, x_ma, y_ma, veta_ini)
            # print('kd = ', kd,' | kt = ', kt, ' | cd = ', cd,' | ct = ', ct)
            
            z_ddot = (-z*(kd+kt) - theta*(kt*b_value-kd*p_value) - zdot*(cd+ct) - thetadot*(ct*b_value-cd*p_value) + kd*desl_rd + kt*desl_rt + cd*Vel_rd + ct*Vel_rt - m_total*9.81) /m_total
            theta_ddot = (-z*(b_value*kt-p_value*kd) - theta*(kt*b_value**2+kd*p_value**2) - zdot*(b_value*ct-p_value*cd) - thetadot*(ct*b_value**2+cd*p_value**2) + kt*b_value*desl_rt - kd*p_value*desl_rd + ct*b_value*Vel_rt - cd*p_value*Vel_rd )/I_total
            # print('z_ddot = ', z_ddot, 'theta_ddot = ', theta_ddot)
               
            return z_ddot, theta_ddot
        
        y = [0,0,0,0]
        dt = 1/600
        z_array = np.zeros_like(t_eval)
        zdot_array = np.zeros_like(t_eval)
        theta_array = np.zeros_like(t_eval)
        thetadot_array = np.zeros_like(t_eval)
        zddot_array = np.zeros_like(t_eval)
        thetaddot_array = np.zeros_like(t_eval)
        acc_ma_array = np.zeros_like(t_eval)
        # count = 0

        for i in range(len(t_eval)):
            z, zdot, theta, thetadot = y
            t = t_eval[i]
            
            z_ddot, theta_ddot = acc_2gdl(t,y)
            acc_ma = z_ddot + theta_ddot*b_value
            # if count >= 60000:
            #     messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            #     count = 0

            zdot += z_ddot*dt
            z += zdot*dt
            thetadot += theta_ddot*dt
            theta += thetadot*dt
            
            z_array[i] = z
            zdot_array[i] = zdot
            theta_array[i] = theta
            thetadot_array[i] = thetadot
            zddot_array[i] = z_ddot
            thetaddot_array[i] = theta_ddot
            acc_ma_array[i] = acc_ma
            y = z, zdot, theta, thetadot

            # count = count + 1
        
        # plt.plot(t_eval, zddot_array, color = 'red')
        # plt.plot(t_eval, thetaddot_array, color = 'blue')
        # plt.xlabel('time [s]')
        # plt.ylabel('aceleração [m/s^2]')
        # plt.show()    
        self.df['acc_z_ma_tr'] = acc_ma_array #aceleração no ponto do amortecedor traseiro [m/s^2]

    def Rainflow_count(self):
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return

        # Get the selected Y-axis columns from the checkboxes
        selected_columns = self.get_selected_columns()
        print('You chose to look for data from: ', selected_columns)
        
        if len(selected_columns) > 1:
            messagebox.showinfo("Atenção!", "Selecione apenas 1 conjunto de dados para contagem Rainflow")
            return
        elif len(selected_columns) == 0:
            messagebox.showinfo("Atenção!", "Selecione qual conjunto de dados deseja aplicar a contagem Rainflow")
            return
        else:

            selected_columns = selected_columns[0]
            #Find the index of peaks and valleys
            peaks, _ = find_peaks(self.df[selected_columns])
            valleys, _ = find_peaks(-self.df[selected_columns])
        
            #Creating dedicated dataframes to store the peaks and valleys values from "Data" associated with their index.
            peaks_df = pd.DataFrame({'index': peaks, selected_columns: self.df.iloc[peaks][selected_columns], 'type': 'peak'})
            valleys_df = pd.DataFrame({'index': valleys, selected_columns: self.df.iloc[valleys][selected_columns], 'type': 'valley'})
        
       
            #Build final data Frame by concatenating and sorting by index the peaks and valleys, then dropping the old index list.
            cycles_peaks = pd.concat([peaks_df, valleys_df]).sort_values(by='index').reset_index(drop=True)
            
            # Assuming cycles_peaks is already created and sorted
            max_range = max(cycles_peaks[selected_columns].max(), abs(cycles_peaks[selected_columns].min()))
            threshold = 0.15 * max_range
            filtered_cycles_peaks = cycles_peaks[abs(cycles_peaks[selected_columns]) > threshold]
            
            cycles_generator = rainflow.extract_cycles(filtered_cycles_peaks[selected_columns])
            cycles = list(cycles_generator)

            self.rainflow_matrix= pd.DataFrame(cycles, columns = ['range', 'mean', 'count', 'i_start', 'i_end'])
            self.rainflow_matrix = self.rainflow_matrix.drop(['i_start','i_end'], axis = 1)
            
            #Creating bin Edges
            rainflow_rangebin_edges = np.linspace(self.rainflow_matrix['range'].min(), self.rainflow_matrix['range'].max(), num=20) 
            rainflow_meanbin_edges = np.linspace(self.rainflow_matrix['mean'].min(), self.rainflow_matrix['mean'].max(), num=20)  
            # Create bins
            self.rainflow_matrix['range_bins'] = pd.cut(self.rainflow_matrix['range'], bins=rainflow_rangebin_edges, include_lowest=True)
            self.rainflow_matrix['mean_bins'] = pd.cut(self.rainflow_matrix['mean'], bins=rainflow_meanbin_edges, include_lowest=True)
            #print(self.rainflow_matrix)
            # Group by bins and calculate mean and count
            grouped_ranflow_matrix = self.rainflow_matrix.groupby(['range_bins','mean_bins']).agg(['mean','count']).reset_index()
            print(grouped_ranflow_matrix)
            
            # Rename columns for clarity
            grouped_ranflow_matrix.columns = ['range_bins', 'mean_bins', 'rangebin_mean','range_count','meanbin_mean','mean_count','count_mean','count_count']
            plt.plot(grouped_ranflow_matrix['rangebin_mean'],grouped_ranflow_matrix['range_count'])
            plt.show()
            # grouped_peaks.columns = ['Bin Range', 'Average Value', 'Count']
            # grouped = grouped_filtered_peaks
            
            # print(grouped)
            # plt.figure(figsize=(10,6))
            # plt.bar(grouped['Average Value'],grouped['Count'], color = 'skyblue', label='Count', alpha=0.7)
            # #plt.plot(grouped['Average Value'], grouped['Count'], color='red', linestyle='-', linewidth=2, markersize=8, label='Distribution curve')
            # plt.title('Count vs. Acc peaks')
            # plt.xlabel('Acceleration')
            # plt.ylabel('Event Count')
            # plt.xticks(rotation=45)
            # plt.grid(True)
            # plt.show()
            
    def Histogram_plot(self):
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return

        # Get the selected Y-axis columns from the checkboxes
        selected_columns = self.get_selected_columns()
        print('You chose to look for data from: ', selected_columns)
        
        if len(selected_columns) > 1:
            messagebox.showinfo("Atenção!", "Selecione apenas 1 conjunto de dados para contagem dos eventos")
            return
        elif len(selected_columns) == 0:
            messagebox.showinfo("Atenção!", "Selecione qual conjunto de dados deseja aplicar a contagem de eventos")
            return
        else:

            selected_columns = selected_columns[0]
            #Find the index of peaks and valleys
            peaks, _ = find_peaks(self.df[selected_columns])
            valleys, _ = find_peaks(-self.df[selected_columns])
        
            #Creating dedicated dataframes to store the peaks and valleys values from "Data" associated with their index.
            peaks_df = pd.DataFrame({'index': peaks, selected_columns: self.df.iloc[peaks][selected_columns], 'type': 'peak'})
            valleys_df = pd.DataFrame({'index': valleys, selected_columns: self.df.iloc[valleys][selected_columns], 'type': 'valley'})
               
            #Build final data Frame by concatenating and sorting by index the peaks and valleys, then dropping the old index list.
            cycles_peaks = pd.concat([peaks_df, valleys_df]).sort_values(by='index').reset_index(drop=True)
            
            #Filtering slow reversals
            max_range = max(cycles_peaks[selected_columns].max(), abs(cycles_peaks[selected_columns].min()))
            threshold = 0.15 * max_range
            filtered_cycles_peaks = cycles_peaks[abs(cycles_peaks[selected_columns]) > threshold]
            
            # Assuming cycles_peaks is already created and sorted
            # Creating bin Edges
            filtered_bin_edges = np.linspace(filtered_cycles_peaks[selected_columns].min(), filtered_cycles_peaks[selected_columns].max(), num=20)  # 10 bins
            bin_edges = np.linspace(cycles_peaks[selected_columns].min(), cycles_peaks[selected_columns].max(), num=20)  # 10 bins
            # Create bins
            filtered_cycles_peaks['bins'] = pd.cut(filtered_cycles_peaks[selected_columns], bins=filtered_bin_edges, include_lowest=True)
            cycles_peaks['bins'] = pd.cut(cycles_peaks[selected_columns], bins=bin_edges, include_lowest=True)
            # Group by bins and calculate mean and count
            grouped_filtered_peaks = filtered_cycles_peaks.groupby('bins')[selected_columns].agg(['mean', 'count']).reset_index()
            grouped_peaks = cycles_peaks.groupby('bins')[selected_columns].agg(['mean', 'count']).reset_index()
            
            # Rename columns for clarity
            grouped_filtered_peaks.columns = ['Bin Range', 'Average Value', 'Count']
            grouped_peaks.columns = ['Bin Range', 'Average Value', 'Count']
            grouped = grouped_filtered_peaks
            
            print(grouped)
            plt.figure(figsize=(10,6))
            plt.bar(grouped['Average Value'],grouped['Count'], color = 'skyblue', label='Count', alpha=0.7)
            plt.plot(grouped['Average Value'], grouped['Count'], color='red', linestyle='-', linewidth=2, markersize=8, label='Distribution curve')
            plt.title('Count vs. Acc peaks')
            plt.xlabel('Acceleration [g]')
            plt.ylabel('Event Count')
            plt.xticks(rotation=45)
            plt.grid(True)
            plt.show()
   
    def PSD(self):

        
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return
        
        # Get the selected X-axis column from the dropdown menu
        x_axis_column = self.x_axis_var.get()
        # Get the selected Y-axis columns from the checkboxes
        selected_columns = self.get_selected_columns()
        
        #Getting num of plots
        num_plots = len(selected_columns)
        
        color_list = ['blue', 'red', 'cyan', 'black', 'green', 'magenta', 'yellow']
        
        if len(selected_columns) > 7:
            messagebox.showinfo("Atenção!", "Selecione no máximo 7 conjunto de dados para contagem dos eventos")
            return
        elif len(selected_columns) == 0:
            messagebox.showinfo("Atenção!", "Selecione qual conjunto de dados deseja aplicar a contagem de eventos")
            return
        else:        
            #Getting num of plots
            num_plots = len(selected_columns)
            
            color_list = ['blue', 'red', 'cyan', 'black', 'green', 'magenta', 'yellow']
            
            # Create a figure and subplots
            plt.close('all')
            width_factor = 0.7
            fig = plt.figure(figsize=((root.winfo_screenwidth()/96)*width_factor, 6))
            

            
            for i in range(num_plots):        

                amplitude = self.df[selected_columns[i]].values

                fs = 600
                f , S = welch(amplitude,fs,nperseg = 1024)
                #S = (S**(1/2))*9.81
                plt.semilogy(f, S, color = color_list[i], label = selected_columns[i])
                #plt.ylim([1e-7, 1e2])
                #plt.xlim([0,100])
                plt.xlabel('frequency [Hz]')
                plt.ylabel('Densidade de Potência')
                plt.legend()
            
            # Clear existing canvas and toolbar
            if self.canvas is not None:
                self.canvas.get_tk_widget().destroy()
            if hasattr(self, "toolbar"):
                self.toolbar.destroy()

            # Create the canvas
            self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
            self.canvas.draw()
            self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
            
            # Add navigation toolbar
            self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
            self.toolbar.update()
            self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
            
    def ISO_2631_forOptim(self):
        #Verificando se alguma leitura de dados foi carregada na        
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return
        
        t = 1/600
        x_axis_column = self.x_axis_var.get()
        def pondera_acc(acc,freq):
            #Coeficientes para um polinomio de 10° Graus ajustado para representar tabela 3 de Freq vs Pesos da norma ISO2631
            Wk = pd.DataFrame()
            # Valores retirados da ISO2631 Tabela  3
            Wk['Frequências'] = [ 0.1, 0.125, 0.16, 0.2, 0.25, 0.315, 0.4, 0.5, 0.63, 0.8,   1, 1.25, 1.6,   2, 2.5, 3.15,   4,    5,  6.3,    8,  10, 12.5,  16,  20,  25, 31.5,  40,  50,  63,  80,  100, 125,  160,  200, 250,  315,  400]
            Wk['Pesos'] =       [31.2,  48.6,   79, 121,  182,   263, 352, 418,  459, 477, 482,  484, 494, 531, 631,  804, 967, 1039, 1054, 1036, 988,  902, 768, 636, 513,  405, 314, 246, 186, 132, 88.7,  54, 28.5, 15.2, 7.9, 3.98, 1.95]
            Wk_freq = Wk['Frequências'].values
            Wk_pesos = (Wk['Pesos'].values)/1000
            cs = CubicSpline(Wk_freq, Wk_pesos)
            peso = cs(freq)

            acc_w = acc*peso

            return acc_w
        
        def freq_weighting(y,t):

            fft_y_complete = scipy.fft.fft(y)
            
            N = fft_y_complete.size
            
            freq_y_complete = scipy.fft.fftfreq(N,t)
            
            fft_y = fft_y_complete[1:N//2]
            freq_y = freq_y_complete[1:N//2]
            fft_y_w = np.copy(fft_y_complete)
            fft_y_w[1:N//2] = pondera_acc(fft_y, freq_y)

            if N%2 == 0:
                                    
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2][::-1])
            else:
                    
                fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2+1][::-1])
                    
            y_w = np.real(scipy.fft.ifft(fft_y_w))        
            return y_w
        acc = self.df['acc_z_ma_tr'].values
        acc_w = freq_weighting(acc,t)
                    
        VDV= (np.sum(np.real(acc_w)**4)/self.df[x_axis_column].max())**(1/4)
        return VDV
 
    def ISO_2631(self):
        
        #Verificando se alguma leitura de dados foi carregada na        
        if self.file_name == 'None':
            messagebox.showinfo("Atenção!", "Carregue dados primeiro!")
            return
        
        def pondera_acc(acc,freq):
            #Coeficientes para um polinomio de 10° Graus ajustado para representar tabela 3 de Freq vs Pesos da norma ISO2631
            Wk = pd.DataFrame()
            # Valores retirados da ISO2631 Tabela  3
            Wk['Frequências'] = [ 0.1, 0.125, 0.16, 0.2, 0.25, 0.315, 0.4, 0.5, 0.63, 0.8,   1, 1.25, 1.6,   2, 2.5, 3.15,   4,    5,  6.3,    8,  10, 12.5,  16,  20,  25, 31.5,  40,  50,  63,  80,  100, 125,  160,  200, 250,  315,  400]
            Wk['Pesos'] =       [31.2,  48.6,   79, 121,  182,   263, 352, 418,  459, 477, 482,  484, 494, 531, 631,  804, 967, 1039, 1054, 1036, 988,  902, 768, 636, 513,  405, 314, 246, 186, 132, 88.7,  54, 28.5, 15.2, 7.9, 3.98, 1.95]
            Wk_freq = Wk['Frequências'].values
            Wk_pesos = (Wk['Pesos'].values)/1000
            cs = CubicSpline(Wk_freq, Wk_pesos)
            peso = cs(freq)

            acc_w = acc*peso

            return acc_w
        
        # Get the selected Y-axis columns from the checkboxes
        selected_columns = self.get_selected_columns()
        x_axis_column = self.x_axis_var.get()

        t = 1/600
        if len(selected_columns) > 1:
            messagebox.showinfo("Atenção!", "Selecione no máximo 1 conjunto de dados para cálculo do VDV")
            return
        # elif len(selected_columns) == 0:
        #     messagebox.showinfo("Atenção!", "Selecione qual conjunto de dados deseja aplicar a contagem de eventos")
        #     return
        else:        
      
        #Coletando leitura de aceleração selecionada pelo usuário e convertendo para [m/s^2]
            acc_medido = (self.df[selected_columns[0]].values)*9.81 
            
            def freq_weighting(y,t):

                fft_y_complete = scipy.fft.fft(y)
            
                N = fft_y_complete.size
            
                freq_y_complete = scipy.fft.fftfreq(N,t)
            
                fft_y = fft_y_complete[1:N//2]
                freq_y = freq_y_complete[1:N//2]
                fft_y_w = np.copy(fft_y_complete)
                fft_y_w[1:N//2] = pondera_acc(fft_y, freq_y)

                if N%2 == 0:
                                    
                    fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2][::-1])
                else:
                    
                    fft_y_w[-N//2+1:] = np.conj(fft_y_w[1:N//2+1][::-1])
                    
                y_w = np.real(scipy.fft.ifft(fft_y_w))        
                return y_w
            acc_w_medido = freq_weighting(acc_medido,t)
            
            VDV_medido = (np.sum(np.real(acc_w_medido)**4)/self.df[x_axis_column].max())**(1/4)
            
            
            if 'acc_z_ma_tr' in self.df.columns :
                acc_teorico = (self.df['acc_z_ma_tr'].values)
                acc_w_teorico = freq_weighting(acc_teorico,t)
                VDV_teorico = (np.sum(np.real(acc_w_teorico)**4)/self.df[x_axis_column].max())**(1/4)       
                messagebox.showinfo("Resultado", f"Vibration Dose Value Medido (VDV) = {VDV_medido: .2f}\nVibration Dose Value Teórico (VDV) = {VDV_teorico: .2f} ") 
            else:    
                messagebox.showinfo("Resultado", f"Vibration Dose Value Medido (VDV) = {VDV_medido: .2f}\nVibration Dose Value Teórico (VDV): Sem dados para cálculo ") 
                      
    def pre_slice_data(self):    

        # Get the selected X-axis column from the dropdown menu
        x_axis_column = self.x_axis_var.get()
        # Get the selected Y-axis columns from the checkboxes
        selected_columns = self.get_selected_columns()
        #selected_columns = selected_columns.append(x_axis_column)
        print('Checkboxes + dropdown', selected_columns, '+', x_axis_column)
        #Criando estruturas para manipular as "fatias" de dados
        
        self.Dados_selecionados = pd.DataFrame()
        self.GPS_Course_Dist  = pd.DataFrame()
        
        for i in range(len(selected_columns)):
            self.Dados_selecionados[selected_columns] = self.df[selected_columns]
            
        self.Dados_selecionados[x_axis_column] = self.df[x_axis_column]
        
        print('Colunas no data frame criado', self.Dados_selecionados.columns)
        # New window for user inputs
        input_window = tk.Toplevel(root)
        input_window.title("Intervalo de tempo para fatiar dados")
        
        # Create a frame for the input fields and button in the new window
        frame = ttk.Frame(input_window, padding="10")
        frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Field 1
        ttk.Label(frame, text="Tempo Inicial:").grid(row=0, column=0, sticky=tk.W)
        entry1 = ttk.Entry(frame, width=20)
        entry1.grid(row=0, column=1)

        # Field 2
        ttk.Label(frame, text="Tempo Final:").grid(row=1, column=0, sticky=tk.W)
        entry2 = ttk.Entry(frame, width=20)
        entry2.grid(row=1, column=1)

        
        # Function to handle the submit action
        def slice_action():

            # Here you can handle the values, e.g., print, store, or use them
            self.t_ini = float(entry1.get())
            self.t_end = float(entry2.get())
            print(f"t_ini: {self.t_ini}, t_end: {self.t_end}")
            Dados_Fatiados = self.Dados_selecionados[(self.Dados_selecionados[x_axis_column] >= self.t_ini) & (self.Dados_selecionados[x_axis_column] <= self.t_end)]
            #GPS_Course_Dist_Fatiado = self.GPS_Course_Dist[(self.GPS_Course_Dist['Time(s)'] >= self.t_ini) & (self.GPS_Course_Dist['Time(s)'] <= self.t_end)]
        
            #Calculando distância percorrida na fatia escolhida
            #Total_Dist_Fatiado = GPS_Course_Dist_Fatiado['Displacement'].sum()
        
            # Gravando fatia em um CSV
            new_file_name = 'dados_de_' + self.file_name + '_t_ini_'+str(self.t_ini)+'_t_out_'+str(self.t_end)+'.csv'
            Dados_Fatiados.to_csv(new_file_name, index = False)
            #print('Distância na fatia de coleta no pavé: ', Total_Dist_Fatiado,'m')
            # Optionally close the input_window after submitting
            
            input_window.destroy()
            
        # Submit Button
        submit_button = ttk.Button(frame, text="Slice", command=slice_action)
        submit_button.grid(row=2, column=0, columnspan=2)
    
    def populate_dropdown(self):
        options = list(self.df.columns)
        self.x_axis_var.set(options[0])
        menu = self.x_axis_dropdown["menu"]
        menu.delete(0, "end")
        for option in options:
            menu.add_command(label=option, command=lambda value=option: self.x_axis_var.set(value))

    def create_checkboxes(self):
        # Clear any existing checkboxes
        for widget in self.checkbox_vars:
            widget.destroy()
        self.checkboxes.clear()
        self.checkbox_vars.clear()

        # Create a frame inside the canvas
        checkboxes_frame = Frame(self.checkbox_canvas)

        self.checkbox_vars = []
        row = 0
        for column in self.df.columns:
            if column.startswith("Time"):
                continue

            var = IntVar()
            checkbox = Checkbutton(checkboxes_frame, text=column, variable=var)
            checkbox.pack(side="top", anchor="w", padx=10, pady=5)

            self.checkbox_vars.append(var)
            self.checkboxes.append(checkbox)

        # Pack the frame inside the canvas
        checkboxes_frame.update_idletasks()  # Ensure checkboxes_frame size is updated
        self.checkbox_canvas.create_window((0, 0), anchor="nw", window=checkboxes_frame,
                                           width=self.checkbox_canvas.winfo_reqwidth())

        # Configure the scrollbar
        self.checkbox_canvas.config(scrollregion=self.checkbox_canvas.bbox("all"))
        self.scrollbar.config(command=self.checkbox_canvas.yview)

    def plot_data(self): 
        # Get the selected X-axis column from the dropdown menu
        x_axis_column = self.x_axis_var.get()

        # Get the selected Y-axis columns from the checkboxes
        selected_columns = self.get_selected_columns()

        if not selected_columns or not x_axis_column:
            return

        # Calculate the number of rows for the subplots
        num_plots = len(selected_columns)

        # Create a figure and subplots
        plt.close('all')
        width_factor = 0.7
        fig, axes = plt.subplots(num_plots, 1, figsize=((root.winfo_screenwidth()/96)*width_factor, 2 * num_plots), squeeze=False, constrained_layout=True)

        # Flatten the 2D array of subplots to a 1D array for easy indexing
        axes = axes.flatten()
        self.axes_subplot = axes
       
        # Plot each selected column on a separate subplot
        for i, (y_axis_column, ax) in enumerate(zip(selected_columns, axes)):
            ax.plot(self.df[x_axis_column], self.df[y_axis_column], label=y_axis_column, color = 'blue')
            ax.set_xlabel(x_axis_column)
            ax.set_ylabel(y_axis_column)
            ax.legend()
            # Connect the zoom_update method to the xlim_changed event
        #print(axes)
        #print('ax no plot_data é igual a: ', ax)
        # Clear existing canvas and toolbar
        if self.canvas is not None:
            self.canvas.get_tk_widget().destroy()
        if hasattr(self, "toolbar"):
            self.toolbar.destroy()

        # Create the canvas
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)

        # Connect the zoom_update method to the draw_event
        self.canvas.mpl_connect('button_press_event', self.on_press)
        
        self.canvas.mpl_connect('button_release_event', self.on_release)
        
        # Add navigation toolbar
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
               
    def on_press(self, event):
        self.click_count = self.click_count + 1
        #print(vars(event))
        # Get the axis (subplot) that triggered the draw event
        current_ax = event.inaxes.get_ylabel()
        #print('plot em que se apertou o botão: ', current_ax)
        axes = self.axes_subplot
        ax = event.inaxes
        #print('ax no qual é aplicado callback: ', ax)
        # Find the index of the current_ax in the list of all subplots
        selected_columns = self.get_selected_columns()
        subplot_index = list(selected_columns).index(current_ax)

        self.index_zoomed_subplot = subplot_index
        # Disconnect the callback from the previous 'ax'
        if hasattr(self, 'previous_ax') and self.previous_ax is not None:
            ax_ = self.previous_ax
            ax_.callbacks.disconnect(ax_.callbacks.callbacks['xlim_changed'].pop(0, None))
            #print('disconnected callbacks for: ', ax_)
            #print('callback:', ax_.callbacks.callbacks)
            
        self.previous_ax = ax
               
        ax.callbacks.connect('xlim_changed', lambda event, ax=ax: self.zoom_update(event, axes,subplot_index))
        
        #print('connected callbacks for: ', ax)

        # Now, you know which subplot (subplot_index) has been zoomed in
        #print(f"Subplot {subplot_index + 1} has been zoomed in.")
        # You can perform specific actions based on the zoomed subplot if needed       
        if event.inaxes is not None:
            # Store the x and y coordinates when a mouse button is pressed
            self.press_x = event.xdata
            self.press_y = event.ydata
            
    def on_release(self, event):
        if hasattr(self, 'press_x') and hasattr(self, 'press_y'):
            # Check if the x and y coordinates are the same (no zoom)
            if event.xdata == self.press_x and event.ydata == self.press_y:
                # Find the corresponding index in the DataFrame based on the coordinates
                selected_index = self.find_index_by_coordinates(self.press_x, self.press_y)

                # Highlight the selected data point on all subplots and update GPS and HTML
                # self.highlight_and_update(selected_index)
                self.click_count = self.click_count - 1
                
            # Reset stored coordinates
            delattr(self, 'press_x')
            delattr(self, 'press_y')
            
    # def highlight_and_update(self, selected_index):
    #     # Clear existing highlighted red dot markers on all subplots
    #     for ax in self.axes_subplot:
    #         for line in ax.lines:
    #             if line.get_markerfacecolor() == 'red':
    #                 line.remove()
    #     # Highlight the selected data point on all subplots
    #     selected_x = self.df.loc[selected_index, self.x_axis_var.get()]

    #     for i, ax in enumerate(self.axes_subplot):
    #         selected_y = self.df.loc[selected_index, self.get_selected_columns()[i]]
    #         ax.plot(selected_x, selected_y, marker='o', markersize=8, color='red', linestyle='None')
        
    #     self.canvas.draw()

    #     # Update the GPS data and HTML file based on the selected data point
    #     self.update_gps_and_html(selected_index)
        
    # #def update_gps_and_html(self, selected_index):
    #     # Extract GPS coordinates based on the selected index
    #     selected_latitude = self.df.loc[selected_index, 'GPS_Latitude']
    #     selected_longitude = self.df.loc[selected_index, 'GPS_Longitude']

    #     # Update the HTML file with the new GPS coordinates
    #     self.plot_gps_data(selected_latitude, selected_longitude)

    #     # You may also perform additional actions based on the selected data point
    #     # ...

    #     print(f"Selected Data Point: Index={selected_index}, Latitude={selected_latitude}, Longitude={selected_longitude}")             
        
    def zoom_update(self, event, axes,subplot_index):
        self.syncro_zooms = self.syncro_zooms + 1
        print('vezes que o zoom foi sincronizado: ', self.syncro_zooms)
        clicks = self.click_count
        print('clicks realizados: ', clicks)
        if clicks >= self.syncro_zooms:
         
            # Disconnect the callbacks before updating to avoid recursion
            axes = self.axes_subplot
            #ax.callbacks.disconnect(ax.callbacks.callbacks['xlim_changed'].pop(subplot_index, None))
            #print('evento passado para zoom_update: ', vars(event))
            # Update the x-axis limits for all subplots based on the first subplot
            print(f"xlim do subplot {subplot_index + 1} de teste: ",axes[subplot_index] , 'é igual a ', axes[subplot_index].get_xlim())
            
            for ax in axes[0:]:
                if ax != axes[subplot_index]:
                    ax.set_xlim(axes[subplot_index].get_xlim())
                    print('xlim do subplot: ',list(axes).index(ax) + 1,' ', ax , 'é igual a ', ax.get_xlim())
            
            self.canvas.draw()
        else:
            self.syncro_zooms = 0
            self.click_count = 0
            
    def find_index_by_coordinates(self, x_coord, y_coord):
        # Find the index of the closest data point based on x_coord and y_coord
        # You may need to adjust this based on your specific data and plot settings
        distances = ((self.df[self.x_axis_var.get()] - x_coord) ** 2 +
                    (self.df[self.get_selected_columns()[self.index_zoomed_subplot]] - y_coord) ** 2)
        selected_index = distances.idxmin()
        return selected_index            
    
    def get_selected_columns(self):
        selected_columns = []

        for var, checkbox in zip(self.checkbox_vars, self.checkboxes):
            if var.get() == 1:  # Check if checkbox is selected
                selected_columns.append(checkbox.cget("text"))

        return selected_columns

    #def plot_gps_data(self, selected_latitude=None, selected_longitude=None):
        # Plot GPS data function in here
        n = len(self.df['GPS_Latitude'])

        if n == 0:
            return None

        centroid_Latidute = self.df['GPS_Latitude'].mean()
        centroid_Longitude = self.df['GPS_Longitude'].mean()
        # define latitude and longitude
        lat, lon = centroid_Latidute, centroid_Longitude

        my_map = folium.Map(location=[lat, lon], zoom_start=18)
        gps_data = self.df[['GPS_Latitude', 'GPS_Longitude']]
        points_data = gps_data.to_dict(orient='records')
        point_count = 0
        for point in points_data:
            if point_count == 0:
                folium.CircleMarker(
                    location=[point["GPS_Latitude"], point["GPS_Longitude"]],
                    radius=5,
                    color='yellow',
                    fill=True,
                    fill_color='yellow',
                    fill_opacity=1
                ).add_to(my_map)
            point_count = point_count + 1
            if point_count == 300:
                point_count = 0
        
        if selected_latitude is not None and selected_longitude is not None:
            # Highlight the selected data point on the map
            folium.Marker(
                location=[selected_latitude, selected_longitude],
                popup="Selected Point",
                icon=folium.Icon(color='red', icon='info-sign')
            ).add_to(my_map)              
                 
        html_file_path = "map_RANDON.html"
        my_map.save(html_file_path)
        # Open the HTML file in a new PyQt window
        self.open_html_in_PyQt(html_file_path)
      
    #def open_html_in_PyQt(self, html_file_path):
        if self.html_viewer is not None:
            self.html_viewer.close()  # Close the existing HTML viewer window

        self.html_viewer = HTMLViewer(html_file_path)
        self.html_viewer.show()

        
if __name__ == "__main__":
    root = Tk()
    # app = QApplication(sys.argv)
    data_viz_app = DataVisualizationApp(root)
    root.mainloop()