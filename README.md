A simple yet useful 2D rectangular plate transient heat conduction simulation using FEM & Implict Finite Difference method.  
This simple MATLAB code is adopted from top88.m by Ole Sigmund et al.  
Currently only support Dirichlet Temperature boundary condition. If you want to simulate under Neumann or even more complex boundary conditions, just feel free to try adjusting this code.  

Here four demonstration simulation examples are presented in videos: (All demos set initial temperature as 273.15K and simulate to 120s)  
Demo1.mp4: 200*200 FEM plate, with 373.15K temperature on the whole left boundary, 273.15K temperature on the whole right boundary;  
Demo2.mp4: 200*200 FEM plate, with a pinned heat flux inflow of 100W/(m^2) at the center;  
Demo3.mp4: 200*200 FEM plate, with heat flux inflow of 10W/(m^2) on the whole left boundary, 273.15K temperature on center of right boundary (1/10 length);  
Demo4.mp4: 200*100 FEM plate, with heat flux inflow of 1W/(m^2) on the whole plate, 273.15K temperature on center of right boundary (1/10 length);  

One more thing--For implementation details, please read FEMdetails.pdf :)