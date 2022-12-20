%piano.m
clear all
%specify the music to be played:
 
%Happy Birthday 
% nmax= 26; %number of notes
% tnote = [0, 0.25, 0.5, 1, 1.5,  2,  3, 3.25, 3.5, 4, 4.5, 5, 6, 6.25, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 9.75, 10, 10.5, 11, 11.5]; %array of onset times of the notes (s)
% dnote = [.2,.2,.4, .4,.4,.9,.2, .2,.4,.4,.4, .9, 0.2, 0.2,.4, .4,.4, .4, .4, 0.4, 0.2, 0.2, .4, .4, .4, 0.9]; %array of durations of the notes (s)
% anote = [ 2, 2, 2, 2,1,1,1.5, 1.5, 1.5,1.5,1, 1, 1, 1, 1, 1,1,1,1, 1, 1,1,1, 1, 1, 1]; %array of relative amplitudes of the notes
% inote = [ 1, 1, 3, 1, 6, 5, 1, 1, 3, 1, 8, 6, 1, 1, 13, 10, 6, 5, 3, 8, 11, 11, 10, 6, 8, 6]+12;
% inote1 = [4, 4, 6, 4, 9, 8, 4, 4, 6, 4, 11, 9];
%array of string indices of the notes
 
%Harry Potter
nmax= 24; %number of notes
tnote = [0, 0.5, 0.5, 1.25, 1.5, 2, 2, 3, 3.5, 3.5, 5, 5, 6.5, 6.5, 7.25, 7.5, 8, 8, 8.5, 10, 10, 11, 11.5, 11.5]; %array of onset times of the notes (s)
dnote = [.4,.7,1.4, .2,.4,.9,1.4, .4, 1.4, 1.4,1.4, 1.4, 0.7, 1.4,.2, .4,.9, 1.4, .4, 1.4, 0.9, 0.4, .9, 1.4]; %array of durations of the notes (s)
anote = [ 2, 2, 0.5, 2, 2, 2,0.5,2,2, 0.5, 2,0.5,2, 0.5, 2, 2,2,0.5,2, 2, 0.5,0.5,2, 0.5]; %array of relative amplitudes of the notes
inote = [ 27, 32, 8, 35, 34, 32, 8, 39, 37, 8, 34, 8, 32, 8, 35, 34, 31, 14, 33, 27, 8, 11, 27, 15];
 
 
 
%initialize string parameters:
L=1;          % length of all strings
J=81;dx=L/(J-1); %number of points/string, space step
flow = 110; %frequency of string with lowest pitch (1/s)
nstrings = 39; %number of strings
for i=1:nstrings
  f(i)=flow*2^((i-1)/12);  % frequency (1/s)
  tau(i)=1.2*(440/f(i));   % decay time (s)
  M(i)=1;                  % mass/length
  T(i)=M(i)*(2*L*f(i))^2;  % tension
  R(i)=(2*M(i)*L^2)/(tau(i)*pi^2); % damping constant 
  %Find the largest stable timestep for string i:
  dtmax(i) = - R(i)/T(i) + sqrt((R(i)/T(i))^2 + dx^2/(T(i)/M(i)));
end
%The timestep of the computation has to be stable for all strings:
dtmaxmin = min(dtmax);
%Now set dt and nskip such that:
%dt<=dtmaxmin, nskip is a positive integer, and dt*nskip = 1/8192.
%Also, make nskip as small as possible, given the above criteria.
nskip = ceil(1/(8192*dtmaxmin));
dt=1/(8192*nskip);
tmax=tnote(nmax)+dnote(nmax); clockmax=ceil(tmax/dt);
%initialize an array that will be used to tell a string 
%when to stop vibrating:
tstop=zeros(nstrings,1);
%initialize arrays to store state of the strings:
H=zeros(nstrings,J);
V=zeros(nstrings,J);
 
%(xh1,xh2)=part of string hit by hammer:
xh1=0.25*L;xh2=0.35*L;
%list of points hit by hammer:
jstrike=ceil(1+xh1/dx):floor(1+xh2/dx);
j=2:(J-1); %list of interior points
%initialize array to store soundwave:
count=0; %initialize count of samples of recorded soundwave
S=zeros(1,ceil(clockmax/nskip)); %array to record soundwave
tsave = zeros(1,ceil(clockmax/nskip)); %array for times of samples
 
n=1 ; %initialize note counter
for clock=1:clockmax
  t=clock*dt;
  while((n<=nmax) && tnote(n)<=t)
    V(inote(n),jstrike)=anote(n); %strike string inote(n)
                                  %with amplitude anote(n)
    tstop(inote(n))=t+dnote(n);   %record future stop time
    n=n+1;                        %increment note counter
  end
  for i=1:nstrings
    if(t > tstop(i))
      H(i,:)=zeros(1,J);
      V(i,:)=zeros(1,J);
    else
      V(i,j)=V(i,j) ...
            +(dt/dx^2)*(T(i)/M(i))*(H(i,j+1)-2*H(i,j)+H(i,j-1)) ...
            +(dt/dx^2)*(R(i)/M(i))*(V(i,j+1)-2*V(i,j)+V(i,j-1));
      H(i,j)=H(i,j)+dt*V(i,j);
    end
  end
  if(mod(clock,nskip)==0)
    count=count+1;
S(count)=sum(H(:,2)); %sample the sound at the present time
    tsave(count)=t;      %record the time of the sample
  end
end
soundsc(S(1:count)) %listen to the sound
plot(tsave(1:count),S(1:count)) %plot the soundwave



