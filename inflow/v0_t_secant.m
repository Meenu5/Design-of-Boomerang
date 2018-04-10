function y = v0_t_secant(u,w)

 flag=0;
 a1=5;b1=15;
 c = (a1*func11(u,w,b1) - b1*func11(u,w,a1))/(func11(u,w,b1) - func11(u,w,a1));
while abs(func11(u,w,c)) > 0.001
    a1 = b1;
    b1 = c;
 c = (a1*func11(u,w,b1) - b1*func11(u,w,a1))/(func11(u,w,b1) - func11(u,w,a1));
    
    flag = flag + 1;

    if(flag == 200)
        break;
    end
    
display(c);
y = c;

end
 
end