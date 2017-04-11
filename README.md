# Refpropsi
This function is a simple wrapper of refprop for Matlab.
Unlike refpropm.m, it uses SI units and is more flexible for the ouput format.
Refpropsi supports only basic properties, please use refpropm for the others.

For more information type 

> doc refpropsi



# Usage
Copy both files refpropsi.m and rp_proto64.m in your Matlab path.

The function looks for refprop library in C:\Program files\Refprop\ or
C:\Program files (x86)\Refprop\ on Windows and /opt/refprop on Linux. 
If it is located elsewhere, please edit the function "find_refprop_path" on 
line 326

Call refpropsi as you would call any other matlab function