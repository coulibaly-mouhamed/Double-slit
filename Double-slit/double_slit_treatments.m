vect_dir= [cos(theta) sin(theta)];
c =-X.^(cos(theta))-Y.^(sin(theta));
pos_mur = -c - 40*sin(theta);
pos_mur = pos_mur./cos(theta);
histogram(real(pos_mur),70);