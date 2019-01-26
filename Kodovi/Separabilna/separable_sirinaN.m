Nf_min = 2;
Nf_max = 16;

sirine = zeros(Nf_max-Nf_min,1);

for Nf=Nf_min:Nf_max
    stanje=zeros(Nf, 2);
    rezultati=zeros(2*Nf+1,1);
    sq3=sqrt(3);

    for i=1:Nf
        stanje(i,1) = rand();
        stanje(i,2) = sqrt(1-stanje(i,1)*stanje(i,1));
    end

    varijacije = pick(-1:1, Nf, 'or');

    for i=1:size(varijacije,1)
        zbir=0;
        p=1;
        for j=1:size(varijacije, 2)
            t=varijacije(i,j);
            zbir=zbir+t;
            a=stanje(j,1);
            b=stanje(j,2);

            if (t==1)
                p=p*(1.0/3 + sq3/6*(a*a-b*b) - 1.0/3*a*b);
            elseif (t==-1) 
                p=p*(1.0/3 - sq3/6*(a*a-b*b) - 1.0/3*a*b);
            else
                p=p*(1.0/3*(1+2*a*b));
            end
        end

        rezultati(zbir+Nf+1,1)=rezultati(zbir+Nf+1,1)+p;
    end

    x=(-Nf:Nf)/Nf;
    %sirine(Nf-Nf_min+1,1)=fwhm(x,rezultati);
    sirine(Nf-Nf_min+1,1)=std(rezultati);
end

X=Nf_min:Nf_max;
g=1./sqrt(X);
plot(X, sirine);
hold on
plot(X, g);