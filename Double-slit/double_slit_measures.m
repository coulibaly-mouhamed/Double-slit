n=50;
x =linspace(1.5,4.5,n);
theta =[];
X=[];
Y=[];
parfor i=1:n
    i
    [theta_r,x_f,y_f]=double_slit(x(i));
    theta =[theta theta_r];
    X = [X x_f];
    Y =[Y y_f];
end

save(['double_slit_measures','.mat'],'theta','X','Y');
