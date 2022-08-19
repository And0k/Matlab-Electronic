function R= R_shmidt_trig(Upower, U1, U2, Ros)
R(2)= Ros*(U2 - U1)/(Upower-U2);
R(1)= II(Ros,R(2))*(Upower - U1)/U1;