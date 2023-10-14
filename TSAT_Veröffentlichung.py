from tkinter import *
from tkinter import filedialog as fd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import base64
import re
import collections
import sqlite3 
from sqlite3 import Error
import time
import datetime
import sys
import threading
from TSAT_logo import *

#Main function tkinter
TA = Tk()
#Fullscreen
TA.state("zoomed")

# Titel, Background and weighted growth of window
TA.title("TSAT")
TA.configure(background="black")
for x in range(11):
    Grid.rowconfigure(TA, x, weight=1)
for x in range(12,15):
    Grid.rowconfigure(TA,x,weight=2)
for x in range(10):
    Grid.columnconfigure(TA,x, weight=1)

#show the logo for 1.5 seconds when starting TSAT
def ShowLogo():
    Logo=TSAT_Logo
    image1 = PhotoImage(data=Logo)
    labelphoto1 = Label (TA, image = image1) 
    labelphoto1.grid(row=0,column=0, columnspan=20, rowspan=20,sticky=N+E+W+S)
    time.sleep(1.5)
    labelphoto1.destroy()
def LogoMultithread():
    runfunction1=threading.Thread(target=ShowLogo)
    runfunction1.start()
LogoMultithread()

#Create new database (pathwindow, function and button)
createdatabasepath= Text(TA, width=80, height=2,wrap=WORD, background="white")
createdatabasepath.grid(row=0,column=1,columnspan=3, sticky=N+E+W+S)
createdatabase=[]
def CreateDatabase():
    del createdatabase[:]
    newdb = fd.asksaveasfilename(initialdir = "/",title = "Select file",defaultextension=".db", filetypes = [("db files", "*.db")])
    conn1 = sqlite3.connect(newdb)
    c=conn1.cursor()
    conn1.commit()
    conn1.close()
    createdatabase.append(newdb)
    createdatabasepath.delete(0.0,END)
    createdatabasepath.insert(END,createdatabase)
Button(TA,text="Create new Database", width=20, command=CreateDatabase).grid(row=0,column=0,sticky=N+E+W+S)

#Connect to an existing database (pathwindow, function and button)
connectdatabasepath= Text(TA, width=80, height=2,wrap=WORD, background="white")
connectdatabasepath.grid(row=1,column=1,columnspan=3, sticky=N+E+W+S)
connectdatabase=[]
def ConnectDatabase():
    del connectdatabase[:]
    newdb2 = fd.askopenfilename(filetypes = [("db files", "*.db")])
    conn2 = sqlite3.connect(newdb2)
    c=conn2.cursor()
    conn2.commit()
    conn2.close()
    connectdatabase.append(newdb2)
    connectdatabasepath.delete(0.0,END)
    connectdatabasepath.insert(END,connectdatabase)
Button(TA,text="Connect to Database", width=20, command=ConnectDatabase).grid(row=1,column=0,sticky=N+E+W+S)

#red TS1 label
TS1label= Label(TA,text="Path TS1:", bg="black", fg="red", font= "none 12 bold").grid(row=2, column=0, sticky=N+E+W+S)
#pathwindow showing the selected TS1 file
TS1pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
TS1pathshown.grid(row=2, column=1, columnspan=3, sticky=N+E+W+S)
#Function to select TS1 file
pathTS1=["none"]
def SelectTS1():
    del pathTS1[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")]) #expand on this in case there is a need for other file formats
    pathTS1.append(file1)
    TS1pathshown.delete(0.0, END)
    TS1pathshown.insert(END,file1)    
Button(TA, text="Select file", width=20, command=SelectTS1).grid(row=2, column=4,sticky=N+E+W+S)
#Variable that checks if TS1 is selected
varTS1 = IntVar()
TS1check= Checkbutton(TA, text= "TS1", variable=varTS1)
TS1check.grid(row=2,column=5,sticky=N+E+W+S)

#Same for TS2
TS2label= Label(TA,text="Path TS2:", bg="black", fg="red", font="none 12 bold").grid(row=3,column=0,sticky=N+E+W+S)
TS2pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
TS2pathshown.grid(row=3, column=1, columnspan=3, sticky=N+E+W+S)
pathTS2=["none"]
def SelectTS2():
    del pathTS2[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathTS2.append(file1)
    TS2pathshown.delete(0.0, END)
    TS2pathshown.insert(END,file1)
Button(TA, text="Select file", width=20, command=SelectTS2).grid(row=3, column=4,sticky=N+E+W+S)
varTS2 = IntVar()
TS2check= Checkbutton(TA, text= "TS2", variable=varTS2)
TS2check.grid(row=3,column=5,sticky=N+E+W+S)

#Same for TS3
TS3label= Label(TA,text="Path TS3:", bg="black", fg="red", font="none 12 bold").grid(row=4,column=0,sticky=N+E+W+S)
TS3pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
TS3pathshown.grid(row=4, column=1, columnspan=3, sticky=N+E+W+S)
pathTS3=["none"]
def SelectTS3():
    del pathTS3[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathTS3.append(file1)
    TS3pathshown.delete(0.0, END)
    TS3pathshown.insert(END,file1)
Button(TA, text="Select file", width=20, command=SelectTS3).grid(row=4, column=4,sticky=N+E+W+S)
varTS3 = IntVar()
TS3check= Checkbutton(TA, text= "TS3", variable=varTS3)
TS3check.grid(row=4,column=5,sticky=N+E+W+S)

#Same for DC2
DC2label= Label(TA,text="Path DC2:", bg="black", fg="red", font="none 12 bold").grid(row=5,column=0,sticky=N+E+W+S)
DC2pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
DC2pathshown.grid(row=5, column=1, columnspan=3, sticky=N+E+W+S)
pathDC2=["none"]
def SelectDC2():
    del pathDC2[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathDC2.append(file1)
    DC2pathshown.delete(0.0, END)
    DC2pathshown.insert(END,file1)    
Button(TA, text="Select file", width=20, command=SelectDC2).grid(row=5, column=4,sticky=N+E+W+S)
varDC2 = IntVar()
DC2check= Checkbutton(TA, text= "DC2", variable=varDC2)
DC2check.grid(row=5,column=5,sticky=N+E+W+S)

#Same for DC3
DC3label= Label(TA,text="Path DC3:", bg="black", fg="red", font="none 12 bold").grid(row=6,column=0,sticky=N+E+W+S)
DC3pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
DC3pathshown.grid(row=6, column=1, columnspan=3, sticky=N+E+W+S)
pathDC3=["none"]
def SelectDC3():
    del pathDC3[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathDC3.append(file1)
    DC3pathshown.delete(0.0, END)
    DC3pathshown.insert(END,file1)
Button(TA, text="Select file", width=20, command=SelectDC3).grid(row=6, column=4,sticky=N+E+W+S)
varDC3 = IntVar()
DC3check= Checkbutton(TA, text= "DC3", variable=varDC3)
DC3check.grid(row=6,column=5,sticky=N+E+W+S)

#Same for ES1
ES1label= Label(TA,text="Path ES1:", bg="black", fg="red", font="none 12 bold").grid(row=7,column=0,sticky=N+E+W+S)
ES1pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
ES1pathshown.grid(row=7, column=1, columnspan=3, sticky=N+E+W+S)
pathES1=["none"]
def SelectES1():
    del pathES1[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathES1.append(file1)
    ES1pathshown.delete(0.0, END)
    ES1pathshown.insert(END,file1)
Button(TA, text="Select file", width=20, command=SelectES1).grid(row=7, column=4,sticky=N+E+W+S)
varES1 = IntVar()
ES1check= Checkbutton(TA, text= "ES1", variable=varES1)
ES1check.grid(row=7,column=5,sticky=N+E+W+S)

#Same for ES2
ES2label= Label(TA,text="Path ES2:", bg="black", fg="red", font="none 12 bold").grid(row=8,column=0,sticky=N+E+W+S)
ES2pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
ES2pathshown.grid(row=8, column=1, columnspan=3, sticky=N+E+W+S)
pathES2=["none"]
def SelectES2():
    del pathES2[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathES2.append(file1)
    ES2pathshown.delete(0.0, END)
    ES2pathshown.insert(END,file1)
Button(TA, text="Select file", width=20, command=SelectES2).grid(row=8, column=4,sticky=N+E+W+S)
varES2 = IntVar()
ES2check= Checkbutton(TA, text= "ES2", variable=varES2)
ES2check.grid(row=8,column=5,sticky=N+E+W+S)

#Same for ES3
ES3label= Label(TA,text="Path ES3:", bg="black", fg="red", font="none 12 bold").grid(row=9,column=0,sticky=N+E+W+S)
ES3pathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
ES3pathshown.grid(row=9, column=1, columnspan=3, sticky=N+E+W+S)
pathES3=["none"]
def SelectES3():
    del pathES3[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathES3.append(file1)
    ES3pathshown.delete(0.0, END)
    ES3pathshown.insert(END,file1)    
Button(TA, text="Select file", width=20, command=SelectES3).grid(row=9, column=4,sticky=N+E+W+S)
varES3 = IntVar()
ES3check= Checkbutton(TA, text= "ES3", variable=varES3)
ES3check.grid(row=9,column=5,sticky=N+E+W+S)

#Same for Library
Liblabel= Label(TA,text="Path Lib:", bg="black", fg="red", font="none 12 bold").grid(row=10,column=0,sticky=N+E+W+S)
Libpathshown = Text(TA, width=80, height=2, wrap=WORD, background="white")
Libpathshown.grid(row=10, column=1, columnspan=3, sticky=N+E+W+S)
pathLib=["none"]
def SelectLib():
    del pathLib[:]
    file1 = fd.askopenfilename(filetypes = [("fastq files", "*.fastq")])
    pathLib.append(file1)
    Libpathshown.delete(0.0, END)
    Libpathshown.insert(END,file1)    
Button(TA, text="Select file", width=20, command=SelectLib).grid(row=10, column=4,sticky=N+E+W+S)
varLib = IntVar()
Libcheck= Checkbutton(TA, text= "Lib", variable=varLib)
Libcheck.grid(row=10,column=5,sticky=N+E+W+S)

#Box Error and Status
currentstatuslabel=Label(TA,text="Current Status:",bg="black",fg="white",font="none 12 bold").grid(row=11, column=0, sticky=N+E+W+S)
currentstatusbox= Text(TA, width=80,height=6,wrap=WORD,background="white")
currentstatusbox.grid(row=11,column=1,columnspan=3,sticky=N+E+W+S)
scrollbar = Scrollbar(TA)
currentstatusbox.config(yscrollcommand=scrollbar.set)
scrollbar.config(command= currentstatusbox.yview)
scrollbar.grid(column=4, row=11, rowspan=1, sticky=N+S+W)

#Choice of Translation direction
translatedirect = []
def translationdirectionselectiondrop(choice): #function that saves the selected direction of Translation
    del translatedirect[:]
    translatedirect.append(choice)

pkvar = StringVar(TA)
choices = {"Forward", "Reverse","Reverse complement","Forward complement", "-----"}
#set the default
pkvar.set ("-----")

popupMenu = OptionMenu(TA, pkvar, *choices, command=translationdirectionselectiondrop) # popupmenu for translation direction selection
popupMenu.grid(row=1, column=4, sticky=N+E+W+S)
Label (TA, text="Choose a translation direction", bg="black" , fg="white", font="none 12 bold") .grid(row=0, column=4, sticky=W)
#Framing amino acids manual entry box and label
Label (TA, text="Select framing aminoacids", bg="black", fg="white", font= "none 12 bold") .grid(row=0, column=5, sticky=N+E+W+S)
FramingAminoAcidsentry = Entry(TA, width=40, bg="white")
FramingAminoAcidsentry.grid(row=1, column=5, sticky=N+E+W+S)

Label (TA, text="or", bg="black", fg="white", font= "none 12 bold") .grid(row=0, column=6, sticky=N+E+W+S)
Label (TA, text="select them from the database", bg="black", fg="white", font= "none 12 bold") .grid(row=0, column=7, sticky=N+E+W+S)
#variable that checks if manual mode is used. If manual mode is used the user needs to provide the framing regions themself and check the box. Framing regions are expected in a form of XXXX(.+?)XXXX
var_manual_framing_region = IntVar()
manual_framing_region_check = Checkbutton(TA, text= "Manual input", variable=var_manual_framing_region)
manual_framing_region_check.grid(row=1,column=6,sticky=N+E+W+S)

dataforprimer= {}
selecteddbname= []
dbframingregion=[]
framingregionchoice={}

def framingdrop(choice1):#saves the translation direction for later use 
    del selecteddbname[:]
    selecteddbname.append(choice1)

def connectprimer():#connects to an existing database with information of framing regions and translation direction or adds new information to that database
    dataforprimer.clear()
    framingregionchoice.clear()
    connectdb2 = fd.askopenfilename(filetypes = [("db files", "*.db")])
    primerdataconnection = sqlite3.connect(connectdb2)
    connprimer=primerdataconnection.cursor()
    connprimer.execute("Select * from all_framing_regions")#get the data from database
    primerdata= connprimer.fetchall()
    for w,x,y,z in primerdata:#save the information in another database for later use
        dataforprimer.update({w:x})
        framingregionchoice.update({w:(x,y,z)})
    primerdataconnection.commit()
    primerdataconnection.close()
    pkvarfr = StringVar(TA)
    popupMenu2 = OptionMenu(TA,pkvarfr,*framingregionchoice,command=framingdrop)
    popupMenu2.grid(row=2, column=7,sticky=N+E+W+S)
    
    def newentry():# new entry into the connected database
        Dat = Toplevel()
        Dat.title("Database new entry")
        Dat.configure(background="black")
        newentryconnect= sqlite3.connect(connectdb2)
        connprimer2=newentryconnect.cursor()
        Label (Dat, text="Enter ID",bg="black", fg="white", font="none 12 bold").grid(row=0, column=2)
        Table_name = Entry(Dat, width=40, bg="white")
        Table_name.grid(row=1, column=2, sticky=W)
        Label (Dat, text="Enter framing regions",bg="black", fg="white", font="none 12 bold").grid(row=0, column=3)
        Table_name1 = Entry(Dat, width=40, bg="white")
        Table_name1.grid(row=1, column=3, sticky=W)
        Label (Dat, text="Enter Author",bg="black", fg="white", font="none 12 bold").grid(row=0, column=4)
        Table_name2 = Entry(Dat, width=40, bg="white")
        Table_name2.grid(row=1, column=4, sticky=W)
        Label (Dat, text="Enter translation direction",bg="black", fg="white", font="none 12 bold").grid(row=0, column=5)
        Table_name3 = Entry(Dat, width=40, bg="white")
        Table_name3.grid(row=1, column=5, sticky=W)
        connprimer2.execute("CREATE TABLE IF NOT EXISTS all_framing_regions(Name TEXT,Sequence TEXT, Added TEXT, Primer TEXT)")
        def insertall():
            ID=Table_name.get()# a name, is displayed in the dropdownmenu.
            framingregions12=Table_name1.get()#the framing regions, expected: XXX(.+?)XXX
            author=Table_name2.get()#who added the particular entry
            Primer12=Table_name.get()#translation direction of new entry
            connprimer2.execute("INSERT INTO all_framing_regions(Name,Sequence,Added,Primer) VALUES(?,?,?,?)",(ID,framingregions12,author,Primer12))
            newentryconnect.commit()
            newentryconnect.close
            Dat.destroy()
        Button (Dat, text="Insert Data", width=20, command=insertall).grid(row=2,column=0,sticky=N+E+W+S)
    
    def refresh(): #refresh the connection to selected database in order to display new entrys
        dataforprimer.clear()
        framingregionchoice.clear()
        primerdataconnection = sqlite3.connect(connectdb2)
        connprimer=primerdataconnection.cursor()
        connprimer.execute("Select * from all_framing_regions")
        primerdata= connprimer.fetchall()
        for w,x,y,z in primerdata:
            dataforprimer.update({w:x})
            framingregionchoice.update({w:(x,y,z)})
        primerdataconnection.commit()
        primerdataconnection.close()
        pkvarfr = StringVar(TA)
        popupMenu2 = OptionMenu(TA,pkvarfr,*framingregionchoice,command=framingdrop)
        popupMenu2.grid(row=2, column=7,sticky=N+E+W+S)
                    
    Button(TA,text="refresh",width=20,command=refresh).grid(row=4,column=7,sticky=N+E+W+S)    
    Button(TA,text="Add new framing regions entry",width=20,command=newentry).grid(row=3,column=7,sticky=N+E+W+S)

Button(TA, text="connect", width=20, command=connectprimer).grid(row=1, column=7,sticky=N+E+W+S)#Button that connects to preexisting database with stored framing regions 

Label (TA, text="User", bg="black", fg="white", font= "none 12 bold") .grid(row=2, column=6, sticky=N+E+W+S)
Userentry = Entry(TA, width=20, bg="white")
Userentry.grid(row=3, column=6, sticky=N+E+W+S)

Label (TA, text="Selection", bg="black", fg="white", font= "none 12 bold") .grid(row=4, column=6, sticky=N+E+W+S)
Selection_name_entry = Entry(TA, width=20, bg="white")
Selection_name_entry.grid(row=5, column=6, sticky=N+E+W+S)

Save_Sequences = []
only_sequences = []
data_only_Sequence = []
data_only_amount = []
normierung = []
NurFasta = []

def startrun():#threading is used in order to display the current status of sequence analysis during running
    runfunction=threading.Thread(target=Sequence_analysis)
    runfunction.start()

def insert_providedData(Round):#feedback to the user which step is currently beeing analyzed
    currentstatusbox.insert(END,"\n" + Round + " checked, trying to analyse provided data")  

def Sequence_analysis():#starts the analysis of provided data, may take a while depending on file size
    Runlist={"TS1":[str(varTS1.get()),pathTS1[0]],"TS2":[str(varTS2.get()),pathTS2[0]],"TS3":[str(varTS3.get()),pathTS3[0]],"DC2":[str(varDC2.get()),pathDC2[0]],"DC3":[str(varDC3.get()),pathDC3[0]],"ES1":[str(varES1.get()),pathES1[0]],"ES2":[str(varES2.get()),pathES2[0]],"ES3":[str(varES3.get()),pathES3[0]],"Lib":[str(varLib.get()),pathLib[0]]}
    currentstatusbox.delete(0.0,END)
    currentstatusbox.insert(END,"Starting Analysis, gathering data")
    for Name,check_path in Runlist.items():
        
        del dbframingregion[:]
        for key,value in dataforprimer.items():
            if key == selecteddbname[0]:
                dbframingregion.append(value)
        
        if int(check_path[0]) == 1:
            t1= threading.Thread(target=insert_providedData, args=(Name,))
            t1.start()
            try:
                file = str(check_path[1]) 
            except IndexError: #exception catcher in case the user forgot to provide a path to the data
                currentstatusbox.insert(END, "\n"+Name+" Checked but no file/wrong file format given!") 
                break
            try:
                if var_manual_framing_region.get() == 1:            # if manual checkbox was selected
                    gesuchte_region = FramingAminoAcidsentry.get()  # data from manual input is used e.g.'ACCTCCACC(.+?)AGAGTGAGA'
                elif var_manual_framing_region.get() == 0:          # else the data from a connected database will be used
                    gesuchte_region = dbframingregion[0]
            except IndexError:
                currentstatusbox.insert(END, "\nERROR no framing regions selected")  #exception in case framing region were not provided
                break
            Rohdaten = []
            currentstatusbox.insert(END,"\nAttempting to open and read "+Name+" file, this may take a while.")
            with open(file,'r',encoding="utf8") as f:
                for line in f:
                    for word in line.split():
                        Rohdaten.append(word)
            currentstatusbox.insert(END,"\n"+Name+" loaded")
            found_Seq = []
            currentstatusbox.insert(END,"\nSearching for Sequences. Please be patient.")
            try:
                for x in Rohdaten:
                    m = re.search(gesuchte_region, x)
                    if m:
                        found_Seq.append(m.group(1))
            except IndexError:
                currentstatusbox.insert(END,"\nSorry no Sequences with these framing regions found!")#exception in case nothing was found
                break
    
            DNA_Seq = []
            zuklein = []
            
            for a in found_Seq:
                if len(a) % 3 == 0:#DNA sequence needs to be divisible by three or no translation is possible (codon triplet)
                    DNA_Seq.append(a)
                elif len(a) %3 != 0:
                    zuklein.append(a)
                    
            currentstatusbox.insert(END,"\nFound " +str(len(DNA_Seq))+ " complete sequences and " +str(len(zuklein))+ " parital sequences") #information for the user about complete sequences and partial sequences
            currentstatusbox.insert(END,"\nSequencs isolated, starting translation")
    
            Protein_Sequence = []
            Forward = ['Forward']
            Reverse = ['Reverse']
            Forwardcomplement = ['Forward complement']
            Reversecomplement = ['Reverse complement']
            empty = ['-----']
                
            if translatedirect == Forward:
                currentstatusbox.insert(END,"\nForward translation selected")
            elif translatedirect == Reverse:
                currentstatusbox.insert(END,"\nReverse translation selected")
            elif translatedirect == Forwardcomplement:
                currentstatusbox.insert(END,"\nForward complement translation selected")
            elif translatedirect == Reversecomplement:
                currentstatusbox.insert(END,"\nReverse complement translation selected")
            #translation occurs after counting to ensure a differentiation based on DNA and not protein sequence                           
            for d in DNA_Seq:
                m = Seq(d)
                if translatedirect == Forward:
                    p = m.translate()
                    Protein_Sequence.append(str(p))
                elif translatedirect == Reverse:
                    rf = m.reverse_complement()
                    f = rf.complement()
                    t = f.translate()
                    Protein_Sequence.append(str(t))
                elif translatedirect == Forwardcomplement:
                    r = m.complement()
                    a = r.translate()
                    Protein_Sequence.append(str(a))
                elif translatedirect == Reversecomplement:
                    b = m.reverse_complement()
                    c = b.translate()
                    Protein_Sequence.append(str(c))
                elif translatedirect == empty:
                    Errorbox.delete(0.0,END)
                    Errorbox.insert(END, "Please select a translation direction!")
                    return
            
            currentstatusbox.insert(END,"\nTranslation finished")            
    
            totallenght = sum(len(i) for i in Protein_Sequence)
    
            Final_Sequences = collections.Counter(Protein_Sequence).most_common()#order translated sequences after most common ones
            del Save_Sequences[:]
            Save_Sequences.append(Final_Sequences)
            del only_sequences[:]
            for key,value in Final_Sequences:
                only_sequences.append(key)
            del data_only_Sequence[:]   
            for key,value in Final_Sequences:
                data_only_Sequence.append(key)
            del data_only_amount[:]   
            for key,value in Final_Sequences:
                data_only_amount.append(value)      
            
            normierung.append(len(DNA_Seq)) #norming the raw data for better comparison between different selection rounds
            
            currentstatusbox.insert(END,"\nFinished analysis, beginning transfer of data to database")
            
            def createandinsert(Input1):#transfer into the selected database
                unix = time.time()
                date = str(datetime.datetime.fromtimestamp(unix).strftime('%Y-%m-%d %H:%M:%S'))
                usertsat=Userentry.get()
                selectionnametsat=Selection_name_entry.get()
                database=sqlite3.connect(str(connectdatabase[0]))
                datacursor= database.cursor()
                datacursor.execute("CREATE TABLE IF NOT EXISTS "+ Input1 + " (User TEXT, Selection TEXT, Sequence TEXT,Amount INTEGER,Länge REAL,Normalized REAL,ppm REAL,Datestamp TEXT)")
                zahl=0
                for d in data_only_Sequence:
                    rawnumber=data_only_amount[zahl]
                    normed = rawnumber/normierung[0]
                    ppm = normed*1000000            #parts per million creation for better comparison between selection rounds
                    datacursor.execute("INSERT INTO " + Input1 + " (User,Selection,Sequence,Amount,Länge,Normalized,ppm,Datestamp) VALUES (?,?,?,?,?,?,?,?)",(usertsat,selectionnametsat,data_only_Sequence[zahl],data_only_amount[zahl],len(data_only_Sequence[zahl]),normed,ppm,date))
                    zahl += 1
                database.commit()
                del normierung[:]            
                database.close()
                
            createandinsert(Name)
            
            currentstatusbox.insert(END,"\nTransfer finished.")
        else:
            continue
   
    currentstatusbox.insert(END,"\nOperation finished.")
#Optional filtering after sequenze analysis is finished
#SQL commands that filter and create new tables according to the selected present selection rounds   
varTS1a = IntVar()#Variables to check which selection rounds are provided
varTS2a = IntVar()
varDC2a = IntVar()
varDC3a = IntVar()
varES1a = IntVar()
varES2a = IntVar()
varES3a = IntVar()
varLiba = IntVar()
varTS3a = IntVar()

TS3checka= Checkbutton(TA, text= "TS3", variable=varTS3a)
ES3checka= Checkbutton(TA, text= "ES3", variable=varES3a)
DC3checka= Checkbutton(TA, text= "DC3", variable=varDC3a)
TS2checka= Checkbutton(TA, text= "TS2", variable=varTS2a)
ES2checka= Checkbutton(TA, text= "ES2", variable=varES2a)
DC2checka= Checkbutton(TA, text= "DC2", variable=varDC2a)
TS1checka= Checkbutton(TA, text= "TS1", variable=varTS1a)
ES1checka= Checkbutton(TA, text= "ES1", variable=varES1a)
Libchecka= Checkbutton(TA, text= "Lib", variable=varLiba)

Label1= Label(TA, text="Please select your tables",bg="black", fg="white", font="none 12 bold") .grid(row=6, column=6, columnspan=2)

TS3checka.grid(row=7,column=6,sticky=N+E+W+S)
TS2checka.grid(row=7,column=7,sticky=N+E+W+S)
TS1checka.grid(row=8,column=6,sticky=N+E+W+S)
DC3checka.grid(row=8,column=7,sticky=N+E+W+S)
DC2checka.grid(row=9,column=6,sticky=N+E+W+S)
ES3checka.grid(row=9,column=7,sticky=N+E+W+S)
ES2checka.grid(row=10,column=6,sticky=N+E+W+S)
ES1checka.grid(row=10,column=7,sticky=N+E+W+S)
Libchecka.grid(row=11,column=6,sticky=N+E+W+S)

länge=[]
names=[]
#possible combinations of selection rounds/ missing combinations are nonsensical
alltables=['TS3', 'ES3', 'DC3', 'TS2', 'ES2', 'DC2', 'TS1', 'ES1', 'Lib']    #x

NoTS3= ['TS2', 'ES2', 'DC2', 'TS1', 'ES1', 'Lib']                            #x
NoDC3= ['TS3', 'ES3', 'TS2', 'ES2', 'DC2', 'TS1', 'ES1', 'Lib']              #x
NoTS2= ['TS3', 'ES3', 'DC3','TS1', 'ES1', 'Lib']                             #x
NoDC2= ['TS3', 'ES3', 'DC3', 'TS2', 'ES2', 'TS1', 'ES1', 'Lib']              #x
NoES2= ['TS3', 'ES3', 'DC3', 'TS2', 'DC2', 'TS1', 'ES1', 'Lib']              #x
NoTS1= ['TS3', 'ES3', 'DC3', 'TS2', 'ES2', 'DC2', 'Lib']                     #x
NoES1= ['TS3', 'ES3', 'DC3', 'TS2', 'ES2', 'DC2', 'TS1', 'Lib']              #x

NoTS3TS2= ['TS1', 'ES1', 'Lib']                                              #x
NoTS3DC2= ['TS2', 'ES2', 'TS1', 'ES1', 'Lib']                                #x
NoTS3TS1= ['TS2', 'ES2', 'DC2', 'Lib']                                       #x
NoTS3ES1= ['TS2', 'ES2', 'DC2', 'TS1', 'Lib']                                #x    
NoDC3TS2= ['TS3', 'ES3', 'TS1', 'ES1', 'Lib']                                #x
NoDC3DC2= ['TS3', 'ES3', 'TS2', 'ES2', 'TS1', 'ES1', 'Lib']                  #x
NoDC3ES2= ['TS3', 'ES3', 'TS2', 'DC2', 'TS1', 'ES1', 'Lib']                  #x
NoDC3TS1= ['TS3', 'ES3', 'TS2', 'ES2', 'DC2', 'Lib']                         #x
NoDC3ES1= ['TS3', 'ES3', 'TS2', 'ES2', 'DC2', 'TS1', 'Lib']                  #x
NoTS2TS1= ['TS3', 'ES3', 'DC3', 'Lib']                                       #x
NoTS2ES1= ['TS3', 'ES3', 'DC3', 'TS1', 'Lib']                                #x
NoDC2ES2= ['TS3', 'ES3', 'DC3', 'TS2', 'TS1', 'ES1', 'Lib']                  #x
NoDC2TS1= ['TS3', 'ES3', 'DC3', 'TS2', 'ES2', 'Lib']                         #x
NoDC2ES1= ['TS3', 'ES3', 'DC3', 'TS2', 'ES2', 'TS1', 'Lib']                  #x
NoES2TS1= ['TS3', 'ES3', 'DC3', 'TS2', 'DC2', 'Lib']                         #x    
NoES2ES1= ['TS3', 'ES3', 'DC3', 'TS2', 'DC2', 'TS1', 'Lib']                  #x

NoTS3DC2TS1= ['TS2', 'ES2', 'Lib']                                           #x
NoTS3DC2ES1= ['TS2', 'ES2', 'TS1', 'Lib']                                    #x
NoDC3TS2TS1= ['TS3', 'ES3', 'Lib']                                           #x
NoDC3DC2ES2= ['TS3', 'ES3', 'TS2', 'TS1', 'ES1', 'Lib']                      #x
NoDC3DC2TS1= ['TS3', 'ES3', 'TS2', 'ES2', 'Lib']                             #x
NoDC3DC2ES1= ['TS3', 'ES3', 'TS2', 'ES2', 'TS1', 'Lib']                      #x   
NoDC2ES2TS1= ['TS3', 'ES3', 'DC3', 'TS2', 'Lib']                             #x
NoDC2ES2ES1= ['TS3', 'ES3', 'DC3', 'TS2', 'TS1', 'Lib']                      #x   

NoDC3DC2ES2TS1=['TS3', 'ES3', 'TS2', 'Lib']                                  #x
NoDC3DC2ES2ES1=['TS3', 'ES3', 'TS2', 'TS1', 'Lib']                           #x
           
def SQliteoperation():
    currentstatusbox.delete(0.0,END)
    currentstatusbox.insert(END, "Starting SQL operation. Please be patient!")
    database=sqlite3.connect(str(connectdatabase[0]))
    #check the state of variables
    dope=(varTS3a.get(),varES3a.get(),varDC3a.get(),varTS2a.get(),varES2a.get(),varDC2a.get(),varTS1a.get(),varES1a.get(),varLiba.get())
    
    for i in dope:
        if i == 1:
            länge.append(i)
        else:
            pass
    
    if varTS3a.get() == 1:
        names.append("TS3")
    else:
        pass
    
    if varES3a.get() == 1:
        names.append("ES3")
    else:
        pass
    
    if varDC3a.get() == 1:
        names.append("DC3")
    else:
        pass

    if varTS2a.get() == 1:
        names.append("TS2")
    else:
        pass

    if varES2a.get() == 1:
        names.append("ES2")
    else:
        pass
    
    if varDC2a.get() == 1:
        names.append("DC2")
    else:
        pass
    
    if varTS1a.get() == 1:
        names.append("TS1")
    else:
        pass

    if varES1a.get() == 1:
        names.append("ES1")
    else:
        pass
    
    if varLiba.get() == 1:
        names.append("Lib")
    else:
        pass
  

    cur= database.cursor()
        
    # first SQL database operation/ all following operations follow the same routine but with different selected selection rounds for comparism
    if names == alltables:#joining the existing tabels onto a new table 
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                    FROM TS3
                    LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                    LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                    LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                    LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                    LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                    LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                    LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence
                    LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
        
        database.commit()
        #ppm values are selected to compare differnt selection rounds
        cur.execute("""SELECT MIN(ppm) 
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                        FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                        FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                        FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                        FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                        FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                        FROM Lib""")
        minLib= cur.fetchone()[0]
        #null and empty values have to be changed
        cur.execute("""UPDATE nonull
                        SET Empty3_ppm = 0
                        WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                        SET Direct3_ppm = 0
                        WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                        SET Target2_ppm = 0
                        WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                        SET Empty2_ppm = 0
                        WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                        SET Direct2_ppm = 0
                        WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                        SET Target1_ppm = 0
                        WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
            
        cur.execute("""UPDATE nonull
                        SET Empty1_ppm = 0
                        WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                        SET Libary = 0
                        WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                        FROM nonull
                        WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        #Emptyfactor is created to boost sequences with 0 occurences in the controls.This pushes the corresponding sequences in the generated empty score.     
        Emptyfactor = 2#if a factor of 2 is to little or too high in can be changed here
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor
        # change values of 0 to the minimum found ppm divided by two for each selection round. This is needed for calculations (0 division)
        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        #only take sequences wiht at least 8 amino acids/ if shorter fragments are wanted change the value three lines down
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        #stop codons are denoted with a *, This configuration is for an amber supressor in which a certain stop codon (UAG) is translated as Glutamine (Q)
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
            
    else:
        pass
    #same for all other possibilites           
    if names == NoTS3:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS2.User, TS2.Selection, TS2.Sequence, (TS2.Amount) Target_Round_2, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS2
                     LEFT JOIN TS1 ON TS1.Sequence = TS2.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS2.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS2.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS2.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS2.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]            
                
        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target2_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target2_ppm/ Empty2_ppm,2) AS Empty_score, ROUND(Target2_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
            
    if names == NoDC3:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoTS2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoDC2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoES2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoTS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
               
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoTS3TS2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS1.User, TS1.Selection, TS1.Sequence, (TS1.Amount) Target_Round_1, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS1
                     LEFT JOIN Lib ON Lib.Sequence = TS1.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS1.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target1_ppm >= Empty1_ppm""")
        database.commit()
        
        Emptyfactor = 2
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target1_ppm/ Empty1_ppm,2) AS Empty_score, ROUND(Target1_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoTS3DC2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS2.User, TS2.Selection, TS2.Sequence, (TS2.Amount) Target_Round_2, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS2
                     LEFT JOIN TS1 ON TS1.Sequence = TS2.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS2.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS2.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS2.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                
        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target2_ppm >= Empty2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target2_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES2 = minES2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target2_ppm/ Empty2_ppm,2) AS Empty_score, ROUND(Target2_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoTS3TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS2.User, TS2.Selection, TS2.Sequence, (TS2.Amount) Target_Round_2, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (Lib.ppm) Libary
                     FROM TS2
                     LEFT JOIN Lib ON Lib.Sequence = TS2.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS2.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS2.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
        
        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target2_ppm/ Empty2_ppm,2) AS Empty_score, ROUND(Target2_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoTS3ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS2.User, TS2.Selection, TS2.Sequence, (TS2.Amount) Target_Round_2, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS2
                     LEFT JOIN TS1 ON TS1.Sequence = TS2.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS2.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS2.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS2.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target2_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor


        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target2_ppm/ Empty2_ppm,2) AS Empty_score, ROUND(Target2_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoDC3TS2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoDC3DC2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass    
    
    if names == NoDC3ES2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                       
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoDC3TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]               
        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoDC3ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoTS2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
    
    else:
        pass
    
    if names == NoTS2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoDC2ES2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                       
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
        
    if names == NoDC2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoDC2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoES2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (DC2.ppm) Direct2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass


    if names == NoES2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (DC2.ppm) Direct2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence
                     LEFT JOIN DC2 ON DC2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC2""")
        minDC2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct2_ppm = 0
                    WHERE Direct2_ppm IS NULL OR Direct2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Direct2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfDC2 = minDC2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct2_ppm =" + str(halfDC2) + " WHERE Direct2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoTS3DC2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS2.User, TS2.Selection, TS2.Sequence, (TS2.Amount) Target_Round_2, (TS2.ppm) Target3_ppm, (ES2.ppm) Empty2_ppm, (Lib.ppm) Libary
                     FROM TS2
                     LEFT JOIN Lib ON Lib.Sequence = TS2.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS2.Sequence""")
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                
        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target2_ppm >= Empty2_ppm AND Target2_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES2 = minES2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target2_ppm/ Empty2_ppm,2) AS Empty_score, ROUND(Target2_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoTS3DC2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS2.User, TS2.Selection, TS2.Sequence, (TS2.Amount) Target_Round_2, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS2
                     LEFT JOIN TS1 ON TS1.Sequence = TS2.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS2.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS2.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]               
        
        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target2_ppm >= Empty2_ppm AND Target2_ppm >= Target1_ppm AND Target2_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES2 = minES2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target2_ppm/ Empty2_ppm,2) AS Empty_score, ROUND(Target2_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    
    if names == NoDC3TS2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass  

    if names == NoDC3DC2ES2:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (TS1.ppm) Target1_ppm, (ES1.ppm) Empty1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN ES1 ON ES1.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES1""")
        minES1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Empty1_ppm = 0
                    WHERE Empty1_ppm IS NULL OR Empty1_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Target1_ppm AND Target1_ppm >= Empty1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfES1 = minES1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty1_ppm =" + str(halfES1) + " WHERE Empty1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoDC3DC2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoDC3DC2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (ES2.ppm) Empty2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN ES2 ON ES2.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM ES2""")
        minES2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Empty2_ppm = 0
                    WHERE Empty2_ppm IS NULL OR Empty2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Empty2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfES2 = minES2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Empty2_ppm =" + str(halfES2) + " WHERE Empty2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoDC2ES2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoDC2ES2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (DC3.ppm) Direct3_ppm, (TS2.ppm) Target2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence
                     LEFT JOIN DC3 ON DC3.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM DC3""")
        minDC3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
                
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Direct3_ppm = 0
                    WHERE Direct3_ppm IS NULL OR Direct3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Direct3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfDC3 = minDC3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Direct3_ppm =" + str(halfDC3) + " WHERE Direct3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass

    if names == NoDC3DC2ES2TS1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]
        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
      
    if names == NoDC3DC2ES2ES1:
        cur.execute("""CREATE TABLE IF NOT EXISTS nonull AS SELECT TS3.User, TS3.Selection, TS3.Sequence, (TS3.Amount) Target_Round_3, (TS3.ppm) Target3_ppm, (ES3.ppm) Empty3_ppm, (TS2.ppm) Target2_ppm, (TS1.ppm) Target1_ppm, (Lib.ppm) Libary
                     FROM TS3
                     LEFT JOIN TS2 ON TS2.Sequence = TS3.Sequence
                     LEFT JOIN TS1 ON TS1.Sequence = TS3.Sequence
                     LEFT JOIN ES3 ON ES3.Sequence = TS3.Sequence
                     LEFT JOIN Lib ON Lib.Sequence = TS3.Sequence""")
    
        database.commit()
        
        cur.execute("""SELECT MIN(ppm)
                    FROM ES3""")
        minES3= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS2""")
        minTS2= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM TS1""")
        minTS1= cur.fetchone()[0]
        cur.execute("""SELECT MIN(ppm)
                    FROM Lib""")
        minLib= cur.fetchone()[0]                        
        
        cur.execute("""UPDATE nonull
                    SET Empty3_ppm = 0
                    WHERE Empty3_ppm IS NULL OR Empty3_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target2_ppm = 0
                    WHERE Target2_ppm IS NULL OR Target2_ppm = ''""")
        database.commit()

        cur.execute("""UPDATE nonull
                    SET Target1_ppm = 0
                    WHERE Target1_ppm IS NULL OR Target1_ppm = ''""")
        database.commit()
        
        cur.execute("""UPDATE nonull
                    SET Libary = 0
                    WHERE Libary IS NULL OR Libary = ''""")
        database.commit()

        cur.execute("""CREATE TABLE IF NOT EXISTS actualValue AS SELECT *
                    FROM nonull
                    WHERE Target3_ppm >= Empty3_ppm AND Target3_ppm >= Target2_ppm AND Target2_ppm >= Target1_ppm AND Target3_ppm >= Libary""")
        database.commit()
        
        Emptyfactor = 2
        halfES3 = minES3/Emptyfactor
        halfTS2 = minTS2/Emptyfactor
        halfTS1 = minTS1/Emptyfactor
        halfLib = minLib/Emptyfactor

        cur.execute("UPDATE actualValue SET Empty3_ppm =" + str(halfES3) + " WHERE Empty3_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target2_ppm =" + str(halfTS2) + " WHERE Target2_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Target1_ppm =" + str(halfTS1) + " WHERE Target1_ppm = 0")
        database.commit()

        cur.execute("UPDATE actualValue SET Libary =" + str(halfLib) + " WHERE Libary = 0")
        database.commit()
        
        cur.execute("""CREATE TABLE IF NOT EXISTS Auswertung AS Select * ,ROUND(Target3_ppm/ Empty3_ppm,2) AS Empty_score, ROUND(Target3_ppm/ Libary,2) AS Enrichment_score
                    from actualValue
                    WHERE length(Sequence) >= 8 
                    order by Empty_score desc""")
        database.commit()
        
        cur.execute("""UPDATE Auswertung
                    SET Sequence = REPLACE(Sequence,'*','Q')""")
        database.commit()
        
    else:
        pass
    #providing the newly created data in a Fasta format for further analysis with other programms (e.g. Hammock)    
    if "TS3" in names:
        cur.execute("select Sequence, Target_Round_3,Empty_score from Auswertung")
        column_names= cur.fetchall()
        count_1=1
        savefasta = fd.asksaveasfilename(initialdir = "/",title = "Select file",defaultextension=".fasta",filetypes = (("All files", "*.*"),))
        fasta_file = open(str(savefasta), "w")
        for d,x,y in column_names:
            fasta_file.write(">"+str(count_1)+"|"+str(x)+"|"+str(y)+"\n"+str(d)+"\n")
            count_1 +=1
        fasta_file.close()
    elif "TS3" not in names and "TS2" in names:
        cur.execute("select Sequence, Target_Round_2,Empty_score from Auswertung")
        column_names2= cur.fetchall()
        count_2=1
        savefasta2 = fd.asksaveasfilename(initialdir = "/",title = "Select file",defaultextension=".fasta",filetypes = (("All files", "*.*"),))
        fasta_file = open(str(savefasta2), "w")
        for d,x,y in column_names2:
            fasta_file.write(">"+str(count_2)+"|"+str(x)+"|"+str(y)+"\n"+str(d)+"\n")
            count_2 +=1
        fasta_file.close()
    elif "TS3" not in names and "TS2" not in names and "TS1" in names:
        cur.execute("select Sequence, Target_Round_1,Empty_score from Auswertung")
        column_names3= cur.fetchall()
        count_3=1
        savefasta3 = fd.asksaveasfilename(initialdir = "/",title = "Select file",defaultextension=".fasta",filetypes = (("All files", "*.*"),))
        fasta_file = open(str(savefasta3), "w")
        for d,x,y in column_names3:
            fasta_file.write(">"+str(count_3)+"|"+str(x)+"|"+str(y)+"\n"+str(d)+"\n")
            count_3 +=1
        fasta_file.close()        
    else:
        currentstatusbox.insert(END, "\nError no target selection found")

    del länge[:]
    del names[:]

def Connect_Database():
    runfunction2=threading.Thread(target=SQliteoperation)
    runfunction2.start()


Button (TA, text= "Sqlite Operation", width=20, command=Connect_Database).grid(row=14,column=6,sticky=N+E+W+S)


def exitbutton():
    TA.destroy()


Button(TA, text="Go", width=20, command=startrun).grid(row=14, column=4,sticky=N+E+W+S)
Button(TA, text="Exit", width=20, command=exitbutton).grid(row=14, column=5,sticky=N+E+W+S)

TA.mainloop()#Mainloop of tkinter, needed for functionality