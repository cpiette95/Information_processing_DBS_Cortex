%% Script for generating OU noise inputs

%NoiseOU_1
th = 0.1;
mu = 120;
sig = 14;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1);
rng(10);                
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_2
th = 0.05;
mu = 50;
sig = 7;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(6);                
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_3
th = 0.05;
mu = 60;
sig = 6;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(1);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_4
th = 0.01;
mu = 100;
sig = 4;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(1);                
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_5
th = 0.01;
mu = 100;
sig = 4;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(17);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_6
th = 0.01;
mu = 200;
sig = 5;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1);
rng(10);                
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_7
th = 0.05;
mu = 80;
sig = 8;
dt = 0.05;
t = 0:dt:500;            
x = zeros(1,length(t)-1); 
rng(12);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);

%NOISEOU_8
th = 0.1;
mu = 100;
sig = 16;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(2);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);


%NOISEOU_9
th = 0.1;
mu = 20;
sig = 2;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(23);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);


%NOISEOU_10
th = 0.1;
mu = 15;
sig = 2;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(11);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);


%NOISEOU_11
th = 0.01;
mu = 20;
sig = 1;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(5);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);


%NOISEOU_12
th = 0.01;
mu = 40;
sig = 2;
dt = 0.05;
t = 0:dt:500;             
x = zeros(1,length(t)-1); 
rng(3);                 
for i = 1:length(t)-1
x(i+1) = x(i)+th*(mu-x(i))*dt+sig*sqrt(dt)*randn;
end
figure;
plot(t,x);