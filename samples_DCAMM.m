rng default
p = makedist('Normal','mu',0.5,'sigma',0.15);
p = truncate(p,0,1);
r = random(p,1,1e4);

histogram(r,200)
save('Data','r')