function [x,y,g,l,z,q] = pauliz_hadamard_2_crypto_on(coin1,type)

%computes cryptocurrency Hadamard transforms in Pauli operator space

if nargin ~= 2

	error('Please enter a coin, transform ordering type - 0 for hadamard, 1 for dyadic, default is sequency');

end

%power of twos table
p_two=[1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608];

%read files
Y=xlsread(coin1);
y=length(Y(:,1));

%generate overnight
O=zeros(y,2);
O(1,1)=0;
O(1,2)=0;

for z=2:y

	O(z,2) = Y(z,1)-Y(z,2);
	O(z,1) = Y(z,2)-Y(z-1,1);

end

%generate intraday
N=zeros(y,2);
N(1,1)=0;
N(1,2)=0;

for z=2:y

	N(z,1) = Y(z,1)-Y(z-1,1);
	N(z,2) = Y(z,2)-Y(z-1,2);

end

%generate hadamard vectors
w=length(N(:,1));
x=length(N(:,2));

%user list sizes
if w >= x

	q = w;

else

	q = x;

end

%find transform size
for h=1:length(p_two)

	if q < p_two(h)

		q0 = p_two(h);

        break;

    end

end

%init matrix
m=zeros(q0,1);
n=zeros(q0,1);
o=zeros(q0,1);
p=zeros(q0,1);

%get delta
m(1:w,1)=sign(O(:,1));
n(1:x,1)=sign(O(:,2));
o(1:w,1)=sign(N(:,1));
p(1:x,1)=sign(N(:,2));

%compute matrices
switch type

	case 0

	z = fwht(m,q0,'hadamard');
	q = fwht(n,q0,'hadamard');
	l = fwht(o,q0,'hadamard');
	g = fwht(p,q0,'hadamard');

	case 1

	z = fwht(m,q0,'dyadic');
	q = fwht(n,q0,'dyadic');
	l = fwht(o,q0,'dyadic');
	g = fwht(p,q0,'dyadic');

	otherwise

	z = fwht(m,q0,'sequency');
	q = fwht(n,q0,'sequency');
	l = fwht(o,q0,'sequency');
	g = fwht(p,q0,'sequency');

end

%plot
%subplot(4,1,1);
%[wave,period,scale,coi,sig95]=wt(q.*z);
%wt(q.*z);
%subplot(4,1,2);
%[wave,period,scale,coi,sig95]=wt(l.*g);
%wt(l.*g);
subplot(4,1,1);
wtc(l,g,'ArrowSize',0);
subplot(4,1,2);
wtc(q,z,'ArrowSize',0);
subplot(4,1,3);
wtc(l,z,'ArrowSize',0);
subplot(4,1,4);
wtc(q,g,'ArrowSize',0);

end