  %% u - (u,v,w) , w = ( p,q,r ) , v0 = guess inflow, Cl, Cd = function of alpha, R = 0.269

function yy = func11(u,w,v0)

R=1;
Cl=0.7; Cd=0.8;

y1=zeros(10000);
j=1;
t=0.01;
for r=0:t:R;
 f=@(x) norm(w_b(u,w,r,x,v0))^2*(Cl*cos(alp(u,w,r,x,v0))+Cd*sin(alp(u,w,r,x,v0)));
 %n=200;
 b=2*pi;a=0;
 h=t; % h=(b-a)/n;
 xi=a:h:b;sum1=0;sum2=0;

  for i=3:2:(numel(xi)-2)
  sum1=sum1 + f(xi(i));
  end
for i=2:2:(numel(xi)-1)
  sum2=sum2 + f(xi(i));
end
y1(j)= y1(j) + h/3*(f(xi(1))+2*sum1+4*sum2+f(xi(end)));
j=j+1;
end

h1=1/j;
for i=3:2:(j-3)
  sum1=sum1 + y1(i);
end
for i=2:2:(j-2)
  sum2=sum2 + y1(i);
end
yy1 = h1/3*(y1(1)+2*sum1+4*sum2+y1(j));



  yy1 = 2*yy1;
  disp(yy1);
  %yy1 = 0.0001*yy1;
   yy1 = 5.22315*0.0001*yy1;
   disp(yy1);
  yy2 = yy1/(sqrt(u(1,1)^2 + (u(3,1)-v0)^2));
     disp(yy2);
  yy = yy2 - v0;
end
