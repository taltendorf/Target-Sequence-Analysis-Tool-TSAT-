Objective of TSAT

TSAT is a computational tool for the rapid analysis of NGS data obtained from Phage display and Mirror Image Phage display selections. It is able to extract desired genetic information, translate the information from DNA into protein and in an optional second step filter the obtained information for sequences with a higher probability to bind targets while excluding sequences that may not be target specific due to unspecific binding, amplification advantages or a bias in the library distribution. 




1. Clone repository

   git clone https://github.com/taltendorf/Target-Sequence-Analysis-Tool-TSAT-.git

2. Install TKinter
   
   Debian/Ubuntu: sudo apt-get install python3-tk

   Fedora: dnf install -y python3-tkinter

   Arch: pacman -Syu tk --noconfirm

   REHL/CentosOS: yum install -< python3-tkinter

   OpenSUSE: zypper in -y pyhton-tk

   Mac: brew install python-tk

3. Install Requirements
   
   pip install -r requirements.txt

4. Run TSAT
   
   python3 TSAT_Veroeffentlichung.py
