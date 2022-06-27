      
theta=1;
h1=1/sqrt(2)*[1;1];
        h2=1/sqrt(2)*[1;exp(-1i*theta)];

cvx_begin

variable x(1,1) complex
variable p(2,1)complex
maximize -x;
subject to
p'*p -1<=0;
x-p'*h1*h1'*p<=0;
cvx_end
%%%%%%%%%%%%
