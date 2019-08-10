%When you experiment with different values, change only values of time
%step, total time and the values of c's and moments of Inertia
h = 0.00001;  %Time step
time = 0:h:0.001;  %The last number in this line is the final time
ux0 = 10;
uy0 = 0;
uz0 = 0;
p0 = 1;
q0 = 1;
r0 = 0;
ux = zeros(size(time));
ux(1) = ux0;
uy = zeros(size(time));
uy(1) = uy0;
uz = zeros(size(time));
uz(1) = uz0;
p = zeros(size(time));
p(1) = p0;
q = zeros(size(time));
q(1) = q0;
r = zeros(size(time));
r(1) = r0;
c1 = 0;c2 = 0;c3 = -9.8;c4 = 0;c5 = 0;c6 = 0;  %Change these c's values according to python code and run.
Ix = 0.001;Iy = 0.000001;Iz = 0.001;Jxy = 0; %Change these I's values according to python code and run.
for i = 1:length(time)-1
    k_ux_1 = c1 - q(i)*uz(i) + r(i)*uy(i);
    k_uy_1 = c2 - r(i)*ux(i) + p(i)*uz(i);
    k_uz_1 = c3 - p(i)*uy(i) + q(i)*ux(i);
    k_p_1 = ((c4 + (Jxy/Iy)*c5 - (Iz - Iy)*q(i)*r(i) - (Jxy/Iy)*(Ix - Iz + Iy - Jxy)*p(i)*r(i))*Iy/(Ix*Iy - Jxy*Jxy));
    k_q_1 = c5/Iy + ((Jxy - Ix + Iz)/Iy)*p(i)*r(i) + (Jxy/Iy)*k_p_1;
    k_r_1 = (c6 - (Iy - Ix)*p(i)*q(i) + Jxy*(p(i)*p(i) - q(i)*q(i)))/Iy;
    
    k_ux_2 = c1 - (q(i)+k_q_1/2)*(uz(i)+k_uz_1/2) + (r(i)+k_r_1/2)*(uy(i)+k_uy_1/2);
    k_uy_2 = c2 - (r(i)+k_r_1/2)*(ux(i)+k_ux_1/2) + (p(i)+k_p_1/2)*(uz(i)+k_uz_1/2);
    k_uz_2 = c3 - (p(i)+k_p_1/2)*(uy(i)+k_uy_1/2) + (q(i)+k_q_1/2)*(ux(i)+k_ux_1/2);
    k_p_2 = ((c4 + (Jxy/Iy)*c5 - (Iz - Iy)*(q(i)+k_q_1/2)*(r(i)+k_r_1/2) - (Jxy/Iy)*(Ix - Iz + Iy - Jxy)*(p(i)+k_p_1/2)*(r(i)+k_r_1/2))*Iy/(Ix*Iy - Jxy*Jxy));
    k_q_2 = c5/Iy + ((Jxy - Ix + Iz)/Iy)*(p(i)+k_p_1/2)*(r(i)+k_r_1/2) + (Jxy/Iy)*k_p_2;
    k_r_2 = (c6 - (Iy - Ix)*(p(i)+k_p_1/2)*(q(i)+k_q_1/2) + Jxy*((p(i)+k_p_1/2)*(p(i)+k_p_1/2) - (q(i)+k_q_1/2)*(q(i)+k_q_1/2)))/Iy;
    
    k_ux_3 = c1 - (q(i)+k_q_2/2)*(uz(i)+k_uz_2/2) + (r(i)+k_r_2/2)*(uy(i)+k_uy_2/2);
    k_uy_3 = c2 - (r(i)+k_r_2/2)*(ux(i)+k_ux_2/2) + (p(i)+k_p_2/2)*(uz(i)+k_uz_2/2);
    k_uz_3 = c3 - (p(i)+k_p_2/2)*(uy(i)+k_uy_2/2) + (q(i)+k_q_2/2)*(ux(i)+k_ux_2/2);
    k_p_3 = ((c4 + (Jxy/Iy)*c5 - (Iz - Iy)*(q(i)+k_q_2/2)*(r(i)+k_r_2/2) - (Jxy/Iy)*(Ix - Iz + Iy - Jxy)*(p(i)+k_p_2/2)*(r(i)+k_r_2/2))*Iy/(Ix*Iy - Jxy*Jxy));
    k_q_3 = c5/Iy + ((Jxy - Ix + Iz)/Iy)*(p(i)+k_p_2/2)*(r(i)+k_r_2/2) + (Jxy/Iy)*k_p_3;
    k_r_3 = (c6 - (Iy - Ix)*(p(i)+k_p_2/2)*(q(i)+k_q_2/2) + Jxy*((p(i)+k_p_2/2)*(p(i)+k_p_2/2) - (q(i)+k_q_2/2)*(q(i)+k_q_2/2)))/Iy;
    
    k_ux_4 = c1 - (q(i)+k_q_3)*(uz(i)+k_uz_3) + (r(i)+k_r_3)*(uy(i)+k_uy_3);
    k_uy_4 = c2 - (r(i)+k_r_3)*(ux(i)+k_ux_3) + (p(i)+k_p_3)*(uz(i)+k_uz_3);
    k_uz_4 = c3 - (p(i)+k_p_3)*(uy(i)+k_uy_3) + (q(i)+k_q_3)*(ux(i)+k_ux_3);
    k_p_4 = ((c4 + (Jxy/Iy)*c5 - (Iz - Iy)*(q(i)+k_q_3)*(r(i)+k_r_3) - (Jxy/Iy)*(Ix - Iz + Iy - Jxy)*(p(i)+k_p_3)*(r(i)+k_r_3))*Iy/(Ix*Iy - Jxy*Jxy));
    k_q_4 = c5/Iy + ((Jxy - Ix + Iz)/Iy)*(p(i)+k_p_3)*(r(i)+k_r_3) + (Jxy/Iy)*k_p_4;
    k_r_4 = (c6 - (Iy - Ix)*(p(i)+k_p_3)*(q(i)+k_q_3) + Jxy*((p(i)+k_p_3)*(p(i)+k_p_3) - (q(i)+k_q_3)*(q(i)+k_q_3)))/Iy;
    
    ux(i+1) = ux(i) + (k_ux_1 + 2*k_ux_2 + 2*k_ux_3 + k_ux_4)*h/6;
    uy(i+1) = uy(i) + (k_uy_1 + 2*k_uy_2 + 2*k_uy_3 + k_uy_4)*h/6;
    uz(i+1) = uz(i) + (k_uz_1 + 2*k_uz_2 + 2*k_uz_3 + k_uz_4)*h/6;
    p(i+1) = p(i) + (k_p_1 + 2*k_p_2 + 2*k_p_3 + k_p_4)*h/6;
    q(i+1) = q(i) + (k_q_1 + 2*k_q_2 + 2*k_q_3 + k_q_4)*h/6;
    r(i+1) = r(i) + (k_r_1 + 2*k_r_2 + 2*k_r_3 + k_r_4)*h/6;
    
    %{
    k1 = 1 - uz(i) - q(i);
    l1 = 2- 2*uz(i) - q(i);
    k2 = 1 - (uz(i) + k1/2) - (q(i)+l1/2);
    l2 = 1 - 2*(uz(i) + k1/2) - (q(i)+l1/2);
    k3 = 1 - (uz(i) + k2/2) - (q(i)+l2/2);
    l3 = 1 - 2*(uz(i) + k2/2) - (q(i)+l2/2);
    k4 = 1 - (uz(i) + k3) - (q(i)+l3);
    l4 = 1 - 2*(uz(i) + k3) - (q(i)+l3);
    uz(i+1) = uz(i) + (k1 + 2*k2 + 2*k3 + k4)*h/6;
    q(i+1) = q(i) + (l1 + 2*l2 + 2*l3 + l4)*h/6;
    %}
end
figure;
plot(time,ux);


    