#!/usr/bin/env python
import sys
import os.path
import tkMessageBox
from itertools import combinations
from collections import Counter
from Tkinter import *
from tkFileDialog import askopenfilename, askdirectory

#####################################################################
#####                                                           #####
#####                 GLOBALS and CONFIGURATION                 #####
#####                                                           #####
#####################################################################
# @FASTA_EXTENSION[list]: list of input sequence extension that the program accepts
# @PROTEIN_GROUPS[dictionary]: dictionary of protein groups derived from similarity matrix

FASTA_EXTENSION = ['.fasta', '.fa', '.fas', '.fsa', '.seq']
PROTEIN_GROUPS = {'PAM250': \
                [('?'), ('M','I','V','L'), ('D','N','H','Q','E'), ('F','I','L'), ('S','P','A'), \
                ('S','A','G'), ('Q','K','N'), ('R','H','Q'), ('S','T'), ('R','K','Q'), ('S','N'), \
                ('F','W'), ('R','W'), ('G','D')], \
                'PAM60': \
                [('?'), ('S','A','T'), ('H','N'), ('S','N'), ('N','D'), ('D','E'), ('H','Q'), \
                ('Q','E'), ('Y','F'), ('R','K'), ('M','I'), ('L','M'), ('I','V')], \
                'BLOSUM40': \
                [('?'), ('S','T'), ('S','A'), ('S','Q','N'), ('H','Y','M'), ('D','E'), ('N','D'), \
                ('H','N'), ('W','Y','F'), ('E','Q','K'), ('K','R'), ('R','Q'), ('L','M','I','V'), \
                ('E','K'), ('A','G'), ('L','F','I'), ('V','T')], \
                'BLOSUM62': \
                [('?'), ('S','T'), ('S','A'), ('S','N'), ('H','Y'), ('D','E'), ('N','D'), ('H','N'), \
                ('W','Y','F'), ('E','Q','K'), ('K','R','Q'),('L','M','I','V')], \
                'BLOSUM80': \
                [('?'), ('Q','R','K'), ('Q','E','K'), ('E','D'), ('D','N'), ('Q','H'), ('Y','H'), \
                ('Y','W'), ('Y','F'), ('S','T'), ('S','A'), ('M','V','L','I')]}
#####################################################################



#####################################################################
#####                                                           #####
#####                       USER INTERFACE                      #####
#####                                                           #####
#####################################################################

######################### FUNCTION GUI_main #########################
# setting the elements to the main window
#####################################################################
def GUI_main():
    # Pre-define global variables
    global open_label
    global output_entry
    global output_checkvar1
    global output_checkvar2
    global output_checkvar3
    global result_seq
    global result_output
    global result_log
    global pgroup_value
    global pgroup_radio1
    global pgroup_radio2
    global pgroup_radio3
    global pgroup_radio4
    global pgroup_radio5
    global pgroup_radio6
    filename = None

    # Initiate TK inter
    root = Tk()
    GUI_init(root)

    #################################################################
    #####                   OPEN FILE section                   #####
    #################################################################
    # ====================== SECTION VARIABLES ======================
    # LOCAL  @open_frame[frame]:
    # GLOBAL @open_label[label]:
    # LOCAL  @open_button[button]:
    # ===============================================================
    open_frame = LabelFrame(text='Input File')
    open_label = Label(open_frame,text='Select a file.', relief=FLAT, width=88, anchor=W)
    open_button = Button(open_frame,text='Open...', width=10, command=OpenFile)

    open_frame.grid(row=0,column=0,columnspan=2,padx=5,pady=5,sticky='WE')
    open_label.grid(row=0,column=0,padx=5,pady=5)
    open_button.grid(row=0,column=1,padx=5,pady=5)


    #################################################################
    #####                  OUTPUT PATH section                  #####
    #################################################################
    # ====================== SECTION VARIABLES ======================
    # LOCAL  @output_frame[frame]:
    # GLOBAL @output_entry[entry]:
    # LOCAL  @output_browse[button]:
    # GLOBAL @output_checkvar1[IntVar]:
    # GLOBAL @output_checkvar2[IntVar]:
    # GLOBAL @output_checkvar3[IntVar]:
    # LOCAL  @output_check_frame[frame]:
    # LOCAL  @output_check1[checkbox]:
    # LOCAL  @output_check2[checkbox]:
    # LOCAL  @output_check3[checkbox]:
    # ===============================================================
    output_frame = LabelFrame(text='Output Path')
    output_entry = Entry(output_frame, width=73)
    output_browse = Button(output_frame,text='Browse...', width=10, command=BrowseFile)

    output_frame.grid(row=1,column=0,columnspan=2,padx=5,pady=5,sticky='WE')
    output_frame.pack_propagate(0)
    output_entry.grid(row=0,column=0,padx=5,pady=5)
    output_browse.grid(row=0,column=1,padx=5,pady=5)

    output_check_frame = Frame(output_frame)
    output_checkvar1 = IntVar()
    output_checkvar2 = IntVar()
    output_checkvar3 = IntVar()
    
    output_check1 = Checkbutton(output_check_frame,text = 'Report File', variable = output_checkvar1, onvalue = 1, offvalue = 0)
    output_check1.select()
    output_check2 = Checkbutton(output_check_frame,text = 'Matrix File', variable = output_checkvar2, onvalue = 1, offvalue = 0)
    output_check2.select()
    output_check3 = Checkbutton(output_check_frame,text = 'Gap Info. File', variable = output_checkvar3, onvalue = 1, offvalue = 0)
    output_check3.select()

    output_check_frame.grid(row=0,column=2,padx=5,pady=5,sticky='E')
    output_check1.grid(row=0,column=0,sticky='W')
    output_check2.grid(row=1,column=0,sticky='W')
    output_check3.grid(row=0,column=1,sticky='W')


    #################################################################
    #####                     RESULT section                    #####
    #################################################################
    # ====================== SECTION VARIABLES ======================
    # LOCAL  @result_frame[frame]:
    # GLOBAL @result_seq[text]:
    # GLOBAL @result_output[text]:
    # GLOBAL @result_log[text]:
    # LOCAL  @result_seq_scroll_y[scrollbar]:
    # LOCAL  @result_seq_scroll_x[scrollbar]:
    # LOCAL  @result_output_scroll_y[scrollbar]:
    # LOCAL  @result_output_scroll_x[scrollbar]:
    # LOCAL  @result_log_scroll[scrollbar]:
    # ===============================================================
    result_frame = Frame()
    result_seq = Text(result_frame,width=15,height=25,state=DISABLED,wrap=NONE)
    result_output = Text(result_frame,width=82,height=25,state=DISABLED,wrap=NONE)
    result_log = Text(result_frame,width=100,height=5,state=DISABLED)
    result_seq_scroll_y = Scrollbar(result_frame,command=result_seq.yview)
    result_seq_scroll_x = Scrollbar(result_frame,command=result_seq.xview,orient=HORIZONTAL)
    result_output_scroll_y = Scrollbar(result_frame,command=result_output.yview)
    result_output_scroll_x = Scrollbar(result_frame,command=result_output.xview,orient=HORIZONTAL)
    result_log_scroll = Scrollbar(result_frame,command=result_log.yview)

    result_seq.config(yscrollcommand=result_seq_scroll_y.set)
    result_seq.config(xscrollcommand=result_seq_scroll_x.set)
    result_output.config(yscrollcommand=result_output_scroll_y.set)
    result_output.config(xscrollcommand=result_output_scroll_x.set)
    result_log.config(yscrollcommand=result_log_scroll.set)

    result_frame.grid(row=2,rowspan=2,column=0,padx=5,pady=5,sticky='W')
    result_seq.grid(row=0,column=0,sticky='W')
    result_seq_scroll_y.grid(row=0,column=1,sticky='NWS')
    result_output.grid(row=0,column=2,sticky='W')
    result_output_scroll_y.grid(row=0,column=3,sticky='NES')
    result_seq_scroll_x.grid(row=1,column=0,sticky='WNE')
    result_output_scroll_x.grid(row=1,column=2,sticky='WNE')
    result_log.grid(row=2,column=0,columnspan=3,sticky='W')
    result_log_scroll.grid(row=2,column=1,columnspan=3,sticky='NES')


    #################################################################
    #####                PROTEIN GROUPING section               #####
    #################################################################
    # ====================== SECTION VARIABLES ======================
    # GLOBAL @pgroup_value[string]:
    # LOCAL  @pgroup_frame[frame]:
    # GLOBAL @pgroup_radio1[radio]:
    # GLOBAL @pgroup_radio2[frame]:
    # GLOBAL @pgroup_radio3[frame]:
    # GLOBAL @pgroup_radio4[frame]:
    # GLOBAL @pgroup_radio5[frame]:
    # GLOBAL @pgroup_radio6[frame]:
    # ===============================================================
    pgroup_value = StringVar()

    pgroup_frame = LabelFrame(text='Protein Grouping')
    pgroup_radio1 = Radiobutton(pgroup_frame,text = 'NONE', variable = pgroup_value, value='NONE')
    pgroup_radio2 = Radiobutton(pgroup_frame,text = 'BLOSUM 40', variable = pgroup_value, value='BLOSUM40')
    pgroup_radio3 = Radiobutton(pgroup_frame,text = 'BLOSUM 62', variable = pgroup_value, value='BLOSUM62')
    pgroup_radio4 = Radiobutton(pgroup_frame,text = 'BLOSUM 80', variable = pgroup_value, value='BLOSUM80')
    pgroup_radio5 = Radiobutton(pgroup_frame,text = 'PAM 60', variable = pgroup_value, value='PAM60')
    pgroup_radio6 = Radiobutton(pgroup_frame,text = 'PAM 250', variable = pgroup_value, value='PAM250')

    pgroup_frame.grid(row=2,column=1,padx=5,pady=5,sticky='NW')
    pgroup_radio1.grid(row=0,column=0,sticky='W')
    pgroup_radio1.select()
    pgroup_radio2.grid(row=1,column=0,sticky='W')
    pgroup_radio3.grid(row=2,column=0,sticky='W')
    pgroup_radio4.grid(row=3,column=0,sticky='W')
    pgroup_radio5.grid(row=4,column=0,sticky='W')
    pgroup_radio6.grid(row=5,column=0,sticky='W')


    #################################################################
    #####                   BOTTOM MENU section                 #####
    #################################################################
    # ====================== SECTION VARIABLES ======================
    # LOCAL  @menu_frame[frame]:
    # LOCAL  @menu_run[button]:
    # LOCAL  @menu_save[button]:
    # LOCAL  @menu_help[button]:
    # LOCAL  @menu_about[button]:
    # LOCAL  @menu_exit[button]:
    # ===============================================================
    menu_frame = Frame()
    menu_run = Button(menu_frame,text='RUN!', width=10, command=main)
    menu_save = Button(menu_frame,text='Save Output', width=10, command=SaveOutput)
    menu_help = Button(menu_frame,text='Help', width=10, command=HelpMenu)
    menu_about = Button(menu_frame,text='About', width=10, command=About)
    menu_exit = Button(menu_frame,text='Exit', width=10, command=ExitProgram)

    menu_frame.grid(row=3,column=1,padx=5,pady=5,sticky='SE')
    menu_run.grid(row=0,column=0,pady=10)
    menu_save.grid(row=1,column=0)
    menu_help.grid(row=2,column=0)
    menu_about.grid(row=3,column=0)
    menu_exit.grid(row=4,column=0)


    # STARTS GUI
    mainloop()

######################### FUNCTION GUI_init #########################
# initiate tkInter GUI main window
# ============================== INPUT ==============================
# @root[GUI]: Main window to contain SiMSA program
#####################################################################
def GUI_init(root):
    # Set window title
    root.title("SiMSA")
    # Set window padding
    root["padx"] = 10
    root["pady"] = 10
    # Lock Resizing
    root.resizable(width=FALSE, height=FALSE)

    # # Set icon photo for the window
    # img = PhotoImage(file='icon.gif')
    # root.tk.call('wm', 'iconphoto', root._w, img)

######################### FUNCTION OpenFile #########################
# get user to select, open, and then read the file
#####################################################################
def OpenFile():
    # ======================= GLOBAL VARIABLES ======================
    # @input_content[string]: holds the content of the input file
    # @seq_type[string]: contain the type of the sequence
    #                       'N' for nucleotide sequence
    #                       'P' for protein sequence
    # ===============================================================
    global input_content
    global seq_type

    # Open the file browser
    filename = askopenfilename()
    # Display the full path of selected file
    open_label.config(text=filename)
    # Clear output name and replace with new one
    output_entry.delete(0,'end')
    output_entry.insert(0,'/'.join(filename.split('/')[0:-1]))

    # If the input is alright
    if check_input():
        # Read the file from the given input file
        input_content = fetch_sequences(open_label.cget('text'))
        # Check for sequence type
        seq_type = check_seq_type(input_content)
        # Enable or Disable protein grouping selection
        # according to the type of sequence
        if seq_type == 'N':
            # Disable protein grouping of the sequence is DNA
            pgroup_radio1.select()
            pgroup_radio1.config(state="disabled")
            pgroup_radio2.config(state="disabled")
            pgroup_radio3.config(state="disabled")
            pgroup_radio4.config(state="disabled")
            pgroup_radio5.config(state="disabled")
            pgroup_radio6.config(state="disabled")
        else:
            # Enable protein grouping of the sequence is protein
            pgroup_radio1.config(state="normal")
            pgroup_radio2.config(state="normal")
            pgroup_radio3.config(state="normal")
            pgroup_radio4.config(state="normal")
            pgroup_radio5.config(state="normal")
            pgroup_radio6.config(state="normal")
    
######################## FUNCTION BrowseFile ########################
# let user browse the directory which the output will be written to
#####################################################################
def BrowseFile():
    # ======================= LOCAL VARIABLES =======================
    # @filepath[string]: holds the full path to the selected directory
    # ===============================================================
    # Let user browse for the directory to write output
    filepath = askdirectory()
    # Clear the old path in the output entry box and write the new one in
    output_entry.delete(0,'end')
    output_entry.insert(0,filepath)

######################## FUNCTION SaveOutput ########################
# save the output file to the file system on the computer
#####################################################################
def SaveOutput():
    if output_entry.get() != '':
        # If the output path is selected
        if result_output.get('1.0', END) not in ['',' ',None,'\n']:
            # If SiMSA has already been run in with this input data

            #########################################################
            #####               write REPORT file               #####
            #########################################################
            if output_checkvar1.get() == 1:
                # ================= LOCAL VARIABLES =================
                # @output_rpt[string]: contain full path with name of
                #                      the report file
                # ===================================================
                output_rpt = output_entry.get() + '/report_rpt.txt'

                # Open file and write its content
                with open(output_rpt,'w') as report_file:
                    # =============== LOCAL VARIABLES ===============
                    # @report_file[file]: file to be written
                    # ===============================================
                    # Write the header section of the output file
                    report_write_header(report_file, seq_info['input_content'], seq_info['seq_type'], open_label.cget('text'))
                    # Write identity statistics
                    report_file.write('\nIdentity:\n')
                    write_stats(result_stats['iden_stat'], report_file)
                    # Write similarity statistics if the input sequence is protein
                    if seq_info['seq_type'] == 'P' and result_stats['sim_stat'] != None and result_matrix['sim_mat'] != None:
                        report_file.write('\nSimilarity:\n')
                        write_stats(result_stats['sim_stat'], report_file)

                # Write message to log frame
                insert_text(result_log, 'Output report saved: %s\n' % output_rpt, True)

            #########################################################
            #####               write MATRIX file               #####
            #########################################################
            if output_checkvar2.get() == 1:
                # ================= LOCAL VARIABLES =================
                # @output_mat[string]: contain full path with name of
                #                      the matrix file
                # ===================================================
                output_mat = output_entry.get() + '/report_mat.txt'

                # Open file and write its content
                with open(output_mat,'w') as matrix_file:
                    # =============== LOCAL VARIABLES ===============
                    # @matrix_file[file]: file to be written
                    # ===============================================
                    # Write the header section of the output file
                    report_write_header(matrix_file, seq_info['input_content'], seq_info['seq_type'], open_label.cget('text'))
                    # Write identity result matrix
                    matrix_file.write('\nIdentity:\n')
                    write_matrix(result_matrix['iden_mat'],matrix_file)
                    # Write similarity result matrix if the sequence is protein
                    if seq_info['seq_type'] == 'P' and result_stats['sim_stat'] != None and result_matrix['sim_mat'] != None:
                        matrix_file.write('\nSimilarity:\n')
                        write_matrix(result_matrix['sim_mat'],matrix_file)
                        matrix_file.write('\n')

                # Write message to log frame
                insert_text(result_log, 'Output matrix saved: %s\n' % output_mat, True)

            #########################################################
            #####             write GAP REPORT file             #####
            #########################################################
            if output_checkvar3.get() == 1:
                # ================= LOCAL VARIABLES =================
                # @output_gap[string]: contain full path with name of
                #                      the gap report file
                # ===================================================
                output_gap = output_entry.get() + '/report_gap.txt'

                # Open file and write its content
                with open(output_gap,'w') as gap_file:
                    # =============== LOCAL VARIABLES ===============
                    # @gap_file[file]: file to be written
                    # ===============================================
                    # Write the header section of the output file
                    report_write_header(gap_file, seq_info['input_content'], seq_info['seq_type'], open_label.cget('text'))
                    # Write gap statistics to the file
                    gap_file.write('\n Pruned Gap Position:\n')
                    write_stats(result_stats['gap_stat'],gap_file)
                    gap_file.write('\n\n')
                    # Write gap matrix to the file
                    write_matrix(result_matrix['gap_mat'],gap_file)

                # Write message to log frame
                insert_text(result_log, 'Gap report saved: %s\n' % output_gap, True)
        else:
            # If SiMSA is not run upon the input file yet, then print the error to log frame
            insert_text(result_log, 'ERROR! The calculation has not been done yet.\n', True)
    else:
        # If the output path is not specify, then print the error to log frame
        insert_text(result_log, 'ERROR! Output destination not specified.\n', True)

####################### FUNCTION ClearResults #######################
# clear all displaying results in the output result frame
#####################################################################
def ClearResults():
    result_seq.config(state=NORMAL)
    result_output.config(state=NORMAL)
    result_seq.delete(1.0,END)
    result_output.delete(1.0,END)
    result_seq.config(state=DISABLED)
    result_output.config(state=DISABLED)

####################### FUNCTION insert_text ########################
# insert text into the frame
# ============================== INPUT ==============================
# @text_object[text]: specify which tkInter text object the message
#                     should be insert into
# @message[string]: string to insert into @text_object
# @islog[boolean]: whether the message is log message or not
#####################################################################
def insert_text(text_object, message, islog = False):
    # Insert text into @text_object
    text_object.config(state=NORMAL)
    text_object.insert(INSERT, message)
    text_object.config(state=DISABLED)
    # If this is log message, go to the buttom of the @text_object
    if islog == True:
        text_object.see(END)

######################### FUNCTION HelpMenu #########################
# open the pop-up help menu of SiMSA software
#####################################################################
def HelpMenu():
    # ======================= LOCAL VARIABLES =======================
    # @help_window[window]: window object of this help menu 
    # @var[StringVar]: StringVar object to hold the message to be display
    # @help_string[string]: contain the string that will be printed to the window
    # @help_content[Message]: message object to hold the help string
    # ===============================================================
    # Initialize the help pop-up window
    help_window = Toplevel(padx=10,pady=10)
    help_window.title("How to use SiMSA")
    help_window.minsize(250, 250)
    help_window.resizable(width=FALSE, height=FALSE)

    # Initialize StringVar object
    var = StringVar()
    # Define the message help menu should contain
    help_string = """   1. Select input file.
    (Input file needs to be in fastA format, and should contain ALIGNED sequences.)

    2. Change output path to where you want to save file. And select which file you would like to save from the checkbox.

    3. If the input sequence is protein, select protein grouping.

    4. Press "RUN", the program will display result on the screen.

    5. If you want to save output file, press "Save Output" button."""
    var.set(help_string)
    # Initialize the message object to contain this help message
    help_content = Message(help_window, textvariable=var)
    help_content.grid(row=0,column=0)
    # Display the window
    help_window.mainloop()

########################### FUNCTION About ##########################
# open the pop-up about menu of SiMSA software
#####################################################################
def About():
    # ======================= LOCAL VARIABLES =======================
    # @about_window[window]: window object of this about menu 
    # @var[StringVar]: StringVar object to hold the message to be display
    # @about_content[Message]: message object to hold the about string
    # ===============================================================
    # Initialize the pop-up window
    about_window = Toplevel(padx=10,pady=10)
    about_window.title("About")
    about_window.minsize(250, 250)
    about_window.resizable(width=FALSE, height=FALSE)

    # Set the string to be displayed
    var = StringVar()
    var.set('Created by Bioinformatics Laboratory, Department of Microbiology, \nFaculty of Science, Mahidol University.\n\n\n For bug reports, questions, or suggestions please contact panawun.p@gmail.com')
    # Initialize the message object to contain the about message
    about_content = Message(about_window, textvariable=var)
    about_content.grid(row=0,column=0)
    # Display the window
    about_window.mainloop()

####################### FUNCTION ExitProgram ########################
# exit the program
#####################################################################
def ExitProgram():
    sys.exit()


#####################################################################
#####                                                           #####
#####                       SiMSA PROGRAM                       #####
#####                                                           #####
#####################################################################
def main():
    # Clear all the past results first
    ClearResults()

    # ======================= MAIN VARIABLES ========================
    # GLOBAL @seq_info[dictionary]: holds sequence information for the input sequence
    #                               contain 2 members:
    #                               1. @input_content[dictionary]: contain sequence name and the sequence
    #                               2. @seq_type[string]: 'N' or 'P' to indicate type of sequence
    # GLOBAL @result_matrix[dictionary]: holds each types of matrix calculated from the input
    #                                    contain 3 members:
    #                                    1. sim_mat
    #                                    2. iden_mat
    #                                    3. gap_mat
    # GLOBAL @result_stats[dictionary]: holds each stats calculated from the input
    #                                   contain 3 members:
    #                                   1. sim_stat
    #                                   2. iden_stat
    #                                   3. gap_stat
    # LOCAL  @similarity[list]: contain the similarity calculated for each pair of sequence
    # LOCAL  @identity[list]: contain the identity calculated for each pair of sequence
    # LOCAL  @gap_count[list]: contain the gap_count calculated for each pair of sequence
    # ===============================================================
    global seq_info
    global result_matrix
    global result_stats
    similarity = []
    identity = []
    gap_count = []

    # Check the input and proceed to calculation
    if check_input() and input_content and seq_type:
        # If input is okay and @input_content and @seq_type is available

        # Get the selected protein grouping information if it is protein sequence
        if seq_type == 'P':
            p_groups = fetch_grouping(str(pgroup_value.get()))
        # Match the input sequences into pairs
        seq_pair = pairing(input_content)

        # Write out the current progress to the log
        insert_text(result_log,'Begin calculation...\n', True)

        # Loop through each pair combination of sequence
        # and calculate similarity and identity
        for each_pair in seq_pair:
            # =================== LOCAL VARIABLES ===================
            # @each_pair[dictionary]: contain the information of this combination
            #                         - @pair[tuple]: contain 2 members, which is the index of the sequence
            #                         - @length[int]: contain the total length of sequence that can be calculate
            # @seq1_i[int]: contain the index of first sequence in the pair
            # @seq2_i[int]: contain the index of another sequence in the pair
            # @seq1[string]: contain the sequence of first sequence in the pair
            # @seq2[string]: contain the sequence of another sequence in the pair
            # =======================================================
            # Get the index of sequence in this pair combination
            seq1_i = each_pair['pair'][0]
            seq2_i = each_pair['pair'][1]
            # Get the protein or nucleotide sequence of this pair combination
            seq1 = input_content[seq1_i]['seq']
            seq2 = input_content[seq2_i]['seq']

            # Prune incalculable position in the sequence out
            prune_index = prune_sequence(seq1, seq2, seq_type)
            # Keep the calculable length
            each_pair['length'] = len(seq1) - len(prune_index)
            # Add list of gap index to @gap_count
            gap_count.append(len(prune_index))
        
            # Calculate pairwise identity of these sequence
            identity.append(parewise_identity(seq1, seq2, prune_index))
            # If the sequence is protein and protein grouping is selected
            # Also calculate pairwise similarity
            if seq_type == 'P' and p_groups != 'NONE':
                similarity.append(pairwise_similarity(seq1, seq2, p_groups, prune_index))

        # Write out the current progress to the log
        insert_text(result_log,'Finished calculation.\n', True)

        # Keep the information of the sequence
        seq_info = {'input_content':input_content, 'seq_type':seq_type}
        # Construct and keep the matrix of identity, similarity, and gap
        result_matrix = {'sim_mat':construct_matrix(range(1, len(input_content)+1), seq_pair, similarity),\
                         'iden_mat':construct_matrix(range(1, len(input_content)+1), seq_pair, identity),\
                         'gap_mat':construct_matrix(range(1, len(input_content)+1), seq_pair, gap_count, True)}
        # Calculate and keep the statistic of identity, similarity, and gap
        result_stats = {'sim_stat':statistic(similarity),\
                        'iden_stat':statistic(identity),\
                        'gap_stat':statistic(gap_count, False, True)}

        # Display the result to the output frame
        display_result(seq_info, result_matrix, result_stats)

####################### FUNCTION check_input ########################
# check input file format
#####################################################################
def check_input():
    # ======================= LOCAL VARIABLES =======================
    # @filename[string]: contain the fullpath and filename of the input file
    # @extension[string]: temporary to hold each file extension in @FASTA_EXTENSION
    # ===============================================================
    # Get the input filename
    filename = open_label.cget('text')
    # Check if filename is specify or not
    if filename != 'Select a file.':
        # If the input file is selected
        # Check the file extension
        if any(extension in filename for extension in FASTA_EXTENSION):
            # If the file extension is okay
            return True
        else:
            # If the selected file is not fasta format file
            # Write the error message to log
            insert_text(result_log,'ERROR!! Input file need to be in fasta format.\n',True)
    else:
        # If no input file is selected yet
        # Write the error message to log
        insert_text(result_log,'ERROR!! Input file not specified.\n',True)
    # Return false if there's any error
    return False

######################## FUNCTION read_fasta ########################
# read sequences from the input fasta file (taken from BioPython)
# ============================== INPUT ==============================
# @input_file_obj[file]: contain information of the input fasta file
# ============================= OUTPUT ==============================
# function return the name of the sequence and the sequence itself
#####################################################################
def read_fasta(input_file_obj):
    # ======================= LOCAL VARIABLES =======================
    # @name[string]: contain sequence name as specify in the input file
    # @sequence[string]: contain amino acid or nucleotide sequence
    # ===============================================================

    # Initialize @name and @sequence
    name, sequence = None, []

    # Loop through each line of the input file and read data
    for line in input_file_obj:
        # ======================= LOCAL VARIABLES =======================
        # @line[string]: a temporary contain each line in the file
        # ===============================================================       

        # Strip out all whitespace characters
        line = line.rstrip()
        # Check whether the line is header line or not
        if line.startswith(">"):
            # If this line is header line
            # Return the sequence in case that name has already been declared
            if name: yield (name, ''.join(sequence))
            # If name is not yet declared, that use this line as the name
            # And initialize the list to contain sequence
            name, sequence = line, []
        else:
            # If the line is not header line - it contain sequence information
            # Append the line to sequence list
            sequence.append(line)
    # Return the last sequence in the file
    if name: yield (name, ''.join(sequence))

##################### FUNCTION fetch_sequences ######################
# get sequence from the given input file name
# ============================== INPUT ==============================
# @input_name[string]: name of the input file
# ============================= OUTPUT ==============================
# @input_content[list]: a list of dictionary with name and sequence information
#                       the dictionary contain following key:
#                       'name': name of the sequence
#                       'seq': protein or DNA sequence
#####################################################################
def fetch_sequences(input_name):
    # ======================= LOCAL VARIABLES =======================
    # @seq_count[int]: use to count number of sequence in the file
    # ===============================================================

    # Initiate important local variables
    input_content = [] # list of dictionary | dictionary contain 'name' and 'seq'
    seq_count = 0

    # Write log
    insert_text(result_log,'Reading file: ' + input_name + '\n', True)

    # Open input fasta file and parse its content
    with open(input_name) as input_file_obj:
        # ===================== LOCAL VARIABLES =====================
        # @input_file_obj[file]: contain information of the input fasta file
        # ===========================================================
        for name, seq in read_fasta(fp):
            # =================== LOCAL VARIABLES ===================
            # @name[string]: a temp for name of each sequence in input file
            # @seq[string]: a sequence for each sequence in input file
            # @input_data[dictionary]: a temporary dict to be added to @input_content
            # =======================================================

            # Construct tmp dictionary and eliminate '>' from sequence name
            input_data = {'name':name.translate(None, '>'), 'seq':seq}
            # Append sequence information to @input_content
            input_content.append(input_data)
            # Update number of sequences
            seq_count += 1

    # If in total there's 0 or 1 sequence, then terminate program
    # Since the similarity and identity cannot be calculated
    if seq_count < 2:
        # Write error message to log
        insert_text(result_log,'INPUT ERROR: 2 or more sequences are required.\n', True)

        # Write log
    insert_text(result_log,'Read completed.\n', True)

    # Then return all the input sequences
    return input_content

###################### FUNCTION check_seq_type ######################
# check the input sequences for its type (protein or DNA) and also
# for its length
# ============================== INPUT ==============================
# @input_content[list]: a list of dictionary with name and sequence information
#                       the dictionary contain following key:
#                       'name': name of the sequence
#                       'seq': protein or DNA sequence
# ============================= OUTPUT ==============================
# function return the type of sequence. 'P' for protein and 'N' for DNA
#####################################################################
def check_seq_type(input_content):
    # ======================= LOCAL VARIABLES =======================
    # @types[list]: contain sequence types of all input sequences
    #               the value in each can be 'P' or 'N'
    # ===============================================================
    # Initiate @types to contain the result of sequence type checking
    types = []
    # Print progress message log
    insert_text(result_log,'Checking sequences.\n', True)

    # Check that every input sequences have the same length
    # Then check the sequence type
    if all(len(i['seq'])==len(input_content[0]['seq']) for i in input_content):
        # If all sequences has the same length
        # Check sequence type of all sequence
        for each_input in input_content:
            # =================== LOCAL VARIABLES ===================
            # @each_input[dictionary]: contain sequence informations
            #                          'name': name of the sequence
            #                          'seq': protein or DNA sequence
            # =======================================================

            # Check the sequence types of sequence in current iteration
            if any(c.upper() in 'EFILPQZ' for c in each_input['seq']):
                # sequence is protein if contain any of 'EFILPQZ'
                types.append('P')
            else:
                # recheck the sequence to ensure it really is DNA
                # ================= LOCAL VARIABLES =================
                # @base_count[list]: contain 2 most common bases in the sequence
                # ===================================================

                # Get 2 most common base from the sequence
                base_count = Counter(each_input['seq'].upper()).most_common(2)

                # Recheck the sequence type using the most common base in @base_count
                if all(x[0] in 'ATCGN-.' for x in base_count):
                    # base ATCG is majority of the sequence
                    types.append('N')
                else:
                    types.append('P')

        # Check if every input sequences is of the same type
        if all(each == types[0] for each in types):
            # If every sequences is the same type
            # Write message to log
            insert_text(result_log,'Sequence check completed.\n', True)
            # Then return the sequence type as single character
            return types[0]
        else:
            # If not all sequence is the same type
            # Write error message to log
            insert_text(result_log,'INPUT ERROR: All inputs must consist of same kind of sequences.\n', True)
    else:
        # If not all sequences have the same length
        # Write error message to log
        insert_text(result_log,'INPUT ERROR: All inputs must have the same length.\n', True)
    # Return false if there's any error
    return False

###################### FUNCTION fetch_grouping ######################
# fetch protein groups from the specified matrix choice
# ============================== INPUT ==============================
# @matrix_choice[string]: matrix choice that the user selected
# ============================= OUTPUT ==============================
# @p_groups[list]: each member contain a tuple of proteins in the same group
#####################################################################
def fetch_grouping(group_choice):
    if group_choice == 'NONE':
        # If user does not wish to calculate the similarity
        return 'NONE'
    else:
        # Else return the list contain protein groups accordingly
        return PROTEIN_GROUPS[group_choice]

######################### FUNCTION pairing ##########################
# generate a list of sequence pair combination
# ============================== INPUT ==============================
# @input_content[list]: a list of dictionary with name and sequence information
# ============================= OUTPUT ==============================
# @pair[list]: list of dictionary contain 'pair' and 'length' where:
#              'pair': list of 2 sequence index
#              'length': length of the sequence after check the 
#####################################################################
def pairing(input_content):
    # ======================= LOCAL VARIABLES =======================
    # @all_combination[list]: list of all combinations with no repeat
    # ===============================================================
    # Initiate @pair list
    pair = []
    # generate combination with no repeat
    all_combination = list(combinations(range(len(input_content)), 2))

    # Loop through each combination
    for each_pair in all_combination:
        # ===================== LOCAL VARIABLES =====================
        # @each_pair[tuple]: a temporary for a combination in @all_combination
        # @tmp_dict[dictionary]: a temporary dictionary to hold data for each combination
        # ===========================================================
        # Initialize pair dictionary with the pair and it's total length
        # Initialize 'length' as 0, the value will be assigned after base pruning
        tmp_dict = {'pair':each_pair, 'length':0} 
        # Append new pair to @pair
        pair.append(tmp_dict)

    return pair

###################### FUNCTION prune_sequence ######################
# prune some ambiguous base positions from the sequence, since we
# cannot calculate similarity and identity for this base positions
# ============================== INPUT ==============================
# @seq1[string]: one of the sequence in the pair
# @seq2[string]: another sequence in the pair to be calculate with @seq1
# @seq_type[string]: sequence type, the value is either 'P' or 'N'
# ============================= OUTPUT ==============================
# @pruned_index[list]: list of sequence index to be disregard from the calculation
#####################################################################
def prune_sequence(seq1, seq2, seq_type):
    # ======================= LOCAL VARIABLES =======================
    # @seq_length[int]: length of the sequence before pruning
    # ===============================================================
    # Get the length of the sequence
    seq_length = len(seq1)
    # Initialize a list to contain positions to prune out
    pruned_index = []
    
    for i in range(len(seq1)):
        # ===================== LOCAL VARIABLES =====================
        # @i[int]: a temporary to keep track of sequence position
        # ===========================================================
        if (seq1[i] in '-.?' and seq2[i] in '-.?') or (seq_type == 'P' and 'X?' in [seq1[i], seq2[i]]) or (seq_type == 'N' and 'N?' in [seq1[i], seq2[i]]):
            # prune if position is (gap with gap, ? with ?) 
            # or (one protein sequence contain X or ?) 
            # or (one DNA sequence contain N or ?)
            pruned_index.append(i)

    return pruned_index

#################### FUNCTION parewise_identity #####################
# calculate pairwise identity between 2 sequences
# ============================== INPUT ==============================
# @seq1[string]: one of the sequence in the pair
# @seq2[string]: another sequence in the pair to be calculate with @seq1
# @prune_index[list]: list of sequence index to ignore in the calculation
# ============================= OUTPUT ==============================
# @percent_identity[float]: identity value of the sequence pair in percentage
#####################################################################
def parewise_identity(seq1, seq2, prune_index):
    # ======================= LOCAL VARIABLES =======================
    # @seq_length[int]: length of the sequence before pruning
    # @matched_position[int]: count of the position with identical base
    # ===============================================================

    # if seq1 and seq2 is exactly the same, then identity equals 100%
    if seq1 == seq2:
        return 100

    # Get the length of sequence and initiate @matched_position
    seq_length = len(seq1)
    matched_position = 0

    # Loop through each position in the sequence and see if it is match
    for i in range(seq_length):
        # ===================== LOCAL VARIABLES =====================
        # @i[int]: a temp to keep track of position in the sequence
        # ===========================================================
        if i not in prune_index and seq1[i] == seq2[i]:
            # If the position is not to be pruned and the base is identical
            # Then add it to the matched position
            matched_position += 1

    # Calculate the identity
    percent_identity = float(matched_position)/float(len(seq1)-len(prune_index))*100.00

    return percent_identity

################### FUNCTION pairwise_similarity ####################
# calculate pairwise identity between 2 sequences
# ============================== INPUT ==============================
# @seq1[string]: one of the sequence in the pair
# @seq2[string]: another sequence in the pair to be calculate with @seq1
# @p_groups[string]: list of protein groups from the selected similarity matrix
# @prune_index[list]: list of sequence index to ignore in the calculation
# ============================= OUTPUT ==============================
# @percent_similarity[float]: similarity value of the sequence pair in percentage
#####################################################################
def pairwise_similarity(seq1, seq2, p_groups, prune_index):
    # ======================= LOCAL VARIABLES =======================
    # @seq_length[int]: length of the sequence before pruning
    # @matched_position[int]: count of the position with identical base
    # ===============================================================

    # if seq1 and seq2 is exactly the same, then identity equals 100%
    if seq1 == seq2:
        return 100

    # Get the length of sequence and initiate @matched_position
    seq_length = len(seq1)
    matched_position = 0

    # Loop through each position in the sequence and see if it is match
    for i in range(seq_length):
        # ===================== LOCAL VARIABLES =====================
        # @i[int]: a temp to keep track of position in the sequence
        # ===========================================================
        if i not in prune_index and (any(set([seq1[i].upper(), seq2[i].upper()]).issubset(group) for group in p_groups) or seq1[i] == seq2[i]):
            # If the position is not to be pruned and (both are in the same protein group, or is exactly the same)
            # Then add it to the matched position
            matched_position += 1

    # Calculate the similarity
    percent_similarity = float(matched_position)/float(len(seq1)-len(prune_index))*100.00

    return percent_similarity

######################## FUNCTION statistic #########################
# calculate statistical values: max, min, median, mode, mean for any given @value_list
# ============================== INPUT ==============================
# @value_list[list]: list of value to get statistics from
# @middle_value[boolean]: whether to calculate median, mode, mean, or not
# @non_percentage[boolean]: whether the value in @value_list is percentage value
# ============================= OUTPUT ==============================
# @return_stats[dictionary]: stored the calculated statistic
#####################################################################
def statistic(value_list, middle_value = True, non_percentage = False):
    # ======================= LOCAL VARIABLES =======================
    # @tmp[list]: temporary variable to hold the list of values
    # @most_common[list]: list of tuples, contain 1 most common value in @value_list
    #                     the tuple contains the value at [0] and number of occurance at [1]
    # ===============================================================

    # Return the function with empty value if there's nothing to calculate from
    if value_list == []:
        return None

    # Initiate @return_stats dictionary
    return_stats = {}

    # Sort input value list for mode calculation
    tmp = [] + value_list
    tmp.sort()

    # Calculate median, mode, and mean if @middle_value is True
    if middle_value:
        # Calculate MEAN
        return_stats['Mean'] = sum(tmp)/len(tmp)
        
        # Calculate MEDIAN
        # Check if the value to calculate have odd or even number of members
        if len(tmp) % 2:
            # If the list of value have odd number of members
            if len(tmp) == 1:
                # If there's only 1 member in the list
                return_stats['Median'] = tmp[0]
            else:
                # If there's 3, 5, 7, ... members in the list
                return_stats['Median'] = tmp[(len(tmp)/2) + 1]
        else:
            # If the list of value have even number of members
            return_stats['Median'] = (tmp[len(tmp)/2] + tmp[len(tmp)/2 + 1])/2.0

        # Get MODE
        # Find the most common member in list
        most_common = Counter(tmp).most_common(1)
        # Check if there is most common member or not
        if most_common[0][1] == 1:
            # If there is no most common member, there's no mode
            return_stats['Mode'] = 'NONE'
        else:
            # If there is most common member, return the value with number of occurances
            return_stats['Mode'] = str(round(most_common[0][0], 2)) + '%  with ' + str(most_common[0][1]) + ' occurances.'

    # Get MAX and MIN
    return_stats['Max'] = max(tmp)
    return_stats['Min'] = min(tmp)

    # If the input value is suppose to be percentage
    if not non_percentage:
        # Iterate through the calculated statistics and format
        for key, value in return_stats.iteritems():
            # =================== LOCAL VARIABLES ===================
            # @key[string]: a temp for the key to each @return_stats
            # @value[string/float]: a temp for value i each @return_stats
            # =======================================================
            if key != 'Mode':
                # Format each value to string in [xx.xx%] format
                return_stats[key] = str(round(value, 2)) + '%'

    return return_stats

####################### FUNCTION init_matrix ########################
# initiate the empty matrix to be written into output file
# ============================== INPUT ==============================
# @header[list]: list of headers to each row and column of the matrix
# ============================= OUTPUT ==============================
# @matrix[list]: a 2 dimensional list, each member contain a row of values
#####################################################################
def init_matrix(header):
    # ======================= LOCAL VARIABLES =======================
    # @total_line[int]: total number of lines this matrix should have
    # ===============================================================

    # Initiate the matrix
    matrix = []
    # Calculate total line include header line, seperator line, and value lines
    total_line = len(header) + 2

    # Initialize the matrix with header and filled value cells with '.'
    for line_count in range(total_line):
        # ===================== LOCAL VARIABLES =====================
        # @line_count[int]: a temp to locate the current row in matrix
        # ===========================================================

        # Construct each line of the matrix
        if line_count == 0:
            # If this is the first line of matrix
            # Then it is the header line
            tmp = ['|'] + header

        elif line_count == 1:
            # If this is the second line of matrix
            # Then it is the seperator line
            tmp = ['='*7 for i in range(total_line-1)]

        else:
            # The rest of lines are the value lines
            # Add header at the left and the rest is values
            tmp = [str(header[line_count-2]) + ' |'] + ['.' for i in range(total_line - 2)] # initialize value as '.'

        # Append the line to @matrix
        matrix.append(tmp)

    return matrix

##################### FUNCTION construct_matrix #####################
# construct the string matrix to be printed as an output
# ============================== INPUT ==============================
# @header[list]: list of headers to each row and column of the matrix
# @seq_pair[list]: list of sequence pairs
# @value[list]: list of values to add to the matrix
# @non_percentage[boolean]: whether the value in @value is percentage value
# ============================= OUTPUT ==============================
# @matrix[list]: a 2 dimensional list, each member contain a row of values
#####################################################################
def construct_matrix(header, seq_pair, value, non_percentage = False):
    # If the value list is empty, then there's nothing to construct
    if value == []:
        return None

    # Initiate the matrix
    matrix = init_matrix(header)

    # Loop through each sequence pairing and add values to matrix
    for each_pair in seq_pair:
        # ===================== LOCAL VARIABLES =====================
        # @each_pair[list]: a temp for each sequence pairing
        # @i1[int]: the greater sequence number label
        # @i2[int]: the lesser sequence number label
        # ===========================================================

        # assign each index value so matrix forms into triangle
        i1 = max(each_pair['pair'][0], each_pair['pair'][1])
        i2 = min(each_pair['pair'][0], each_pair['pair'][1])

        # Add the value to matrix position
        if not non_percentage:
            # If the value should be the percentage, add '%' and format
            matrix[i1+2][i2+1] = '%.2f' % round(value[seq_pair.index(each_pair)], 2) + '%'
        else:
            # If the value should not be the percentage, add the value as it is
            matrix[i1+2][i2+1] = '%d' % round(value[seq_pair.index(each_pair)], 2)

    return matrix

###################### FUNCTION display_result ######################
# display the result from calculation in the window
# ============================== INPUT ==============================
# @seq_info[dictionary]: contain all information of the sequence information
# @result_matrix[dictionary]: contain matrix of all calculation
#                               key - 
#                               value - the matrix 
# @result_stats[dictionary]: contain the statistic of all calculation
#####################################################################
def display_result(seq_info, result_matrix, result_stats):
    # ======================= LOCAL VARIABLES =======================
    # @seq_count[int]: use to count the number of sequence for number label
    # ===============================================================
    insert_text(result_seq,'Sequences:\n')

    # Write out each sequence name with number label
    seq_count = 1
    for each_seq in seq_info['input_content']:
        # ===================== LOCAL VARIABLES =====================
        # @each_seq[string]: the name of each sequence
        # ===========================================================
        insert_text(result_seq,'\n%d. %s' % (seq_count,each_seq['name']))
        seq_count += 1
    
    # Write out identity statistic and matrix
    insert_text(result_output, 'Identity:\n')
    write_stats(result_stats['iden_stat'])
    write_matrix(result_matrix['iden_mat'])
    # If the statistic and matrix for similarity exists, also write those out
    if result_stats['sim_stat'] != None and result_matrix['sim_mat'] != None:
        insert_text(result_output, '\n\nSimilarity:\n')
        write_stats(result_stats['sim_stat'])
        write_matrix(result_matrix['sim_mat'])
    # Write out gap statistic and matrix
    insert_text(result_output, '\n\nSequence Gap Information:\n')
    write_stats(result_stats['gap_stat'])
    write_matrix(result_matrix['gap_mat'])

################### FUNCTION report_write_header ####################
# write the header section to the specified output file
# ============================== INPUT ==============================
# @file_data[file]: a file object of the file to write this header to
# @input_content[list]: a list of dictionary with name and sequence information
# @seq_type[string]: sequence type, the value is either 'P' or 'N'
# @input_name[string]: input file name
# ============================= OUTPUT ==============================
# this function write header content to the specified file object
#####################################################################
def report_write_header(file_data, input_content, seq_type, input_name):
    # Write the input file and its path
    file_data.write('Input File: ' + str(input_name) + '\n')

    # Write the type of this sequence
    file_data.write('Sequence Type: ')
    if seq_type == 'P':
        file_data.write('Protein\n')
    else:
        file_data.write('Nucleotide\n')

    # Write names of the input sequences
    file_data.write('\nSequences:\n')
    for each_seq in input_content:
        # ===================== LOCAL VARIABLES =====================
        # @each_seq[dictionary]: a temp contain informations of each sequence
        # ===========================================================

        # Write number labels and sequence name to the file
        file_data.write('    ' + str(input_content.index(each_seq)+1) + '. ' + each_seq['name'] + '\n')

######################## FUNCTION write_stats #######################
# write out the each calculated statistic values to file or output frame
# ============================== INPUT ==============================
# @stats[dictionary]: contain statistical values for each type of stats
#                       key - contain string indicating which tyoe of stats this one is
#                       value - the calculated value of this type of stat
# @file_data[file]: file object to write the information to
#####################################################################
def write_stats(stats, file_data = None):
    # Loop through each stats and print out
    for key, value in stats.iteritems():
        # ===================== LOCAL VARIABLES =====================
        # @key[string]: a temp contain which stats this is
        # @value[float/string]: a temp contain the value of this stats
        # ===========================================================
        if file_data == None:
            # If the file object to write to is not specify
            # Then write the info to the output frame
            insert_text(result_output, '\t%s: %s\n' % (key, value))
        else:
            # If file object is specify
            # Then write the info to the file
            file_data.write('    ' + key + ': ' + str(value) + '\n')

####################### FUNCTION write_matrix #######################
# write out each matrix to the file or to output frame
# ============================== INPUT ==============================
# @matrix[list]: a list of list contain rows and columns of the matrix
# @file_data[file]: file object to write the information to
#####################################################################
def write_matrix(matrix, file_data = None):
    # Loop through the matrix and print out
    for row in matrix:
        # ===================== LOCAL VARIABLES =====================
        # @row[list]: a temp contain information on a row of matrix
        # ===========================================================
        if file_data == None:
            # If the file object to write to is not specify
            insert_text(result_output, '\n')
        else:
            # If file object is specify
            file_data.write('\n')

        for col in row:
            # =================== LOCAL VARIABLES ===================
            # @col[string]: a temp contain each column in the @row
            # =======================================================
            if file_data == None:
                # If the file object to write to is not specify
                # Then write the info to the output frame
                insert_text(result_output,'{:>8s}'.format(str(col)) + ' ')
            else:
                # If file object is specify
                # Then write the info to the file
                file_data.write('{:>8s}'.format(str(col)) + '\t')

#####################################################################
#####                          GLOBAL                           #####
#####################################################################
if __name__ == "__main__":
    GUI_main()