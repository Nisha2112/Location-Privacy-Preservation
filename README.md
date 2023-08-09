# Location-Privacy-Preservation
<h2>General Info</h2>
<p>This project was developed to demonstrate the network model and location preserving routing schemes discussed in the research paper titled “Context-Oriented Location Privacy Preserving Routing Schemes for Wireless Sensor Network”.
</p>
<p>Description of the software functionalities can be found in the paper. A link to the code ocean capsule can also be found in the paper for reference. </p>
<h2>Technology</h2>
<p>The project is created in MATLAB  in version R2018b. MATLAB (an abbreviation of "MATrix LABoratory”) is a proprietary multi-paradigm programming language and numeric computing environment developed by MathWorks. MATLAB allows matrix manipulations, plotting of functions and data, implementation of algorithms, creation of user interfaces, and interfacing with programs written in other languages.</p>
<h1>Installation</h1>
<p>The MATLAB installation steps can be followed from the provided link https://in.mathworks.com/help/install/ug/install-using-a-file-installation-key.html 
</p>
<p>
Follow these steps to install on Linux:
1.	At the system terminal, unzip the matlab_R2023a_glnxa64.zip installer archive to the matlab_R2023a_glnxa64 directory by entering:
unzip matlab_R2023a_glnxa64.zip -d matlab_R2023a_glnxa64
2.	Navigate to the matlab_R2023a_glnxa64 directory.
•	To launch the installer and install to a directory where you have write permissions, execute:
./install
When prompted by the installer, specify the folder for installation.
•	To launch the installer as root, execute:
sudo ./install
If the installer fails to launch as root, it might not have access to the display. Try this workaround:
xhost +SI:localuser:root
sudo -H ./install
xhost -SI:localuser:root
This allows the root user to access the running X server, launches the installer, and then removes the root user from accessing the X server.
</p>

