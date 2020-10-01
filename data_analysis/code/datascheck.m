for i=1:20
y=datasample(x,1000);
std(y)
mean(y)
figure('vis','on')
hist(y)
pause(2)
end