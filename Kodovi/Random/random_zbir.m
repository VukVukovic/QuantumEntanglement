Nf=6;

o0 = [1 0;
      0 1];
ox = [0 1;
      1 0];
oy = [0 complex(0,-1);
      complex(0,1) 0];
oz = [1 0;
      0 -1];

m0 = [1 0 0];
m1 = [-1/2 0 sqrt(3)/2];
m_1 = [-1/2 0 -sqrt(3)/2];

E0 = 1/3*(o0 + m0(1)*ox + m0(2)*oy + m0(3)*oz);
E1 = 1/3*(o0 + m1(1)*ox + m1(2)*oy + m1(3)*oz);
E_1 = 1/3*(o0 + m_1(1)*ox + m_1(2)*oy + m_1(3)*oz);

varijacije = pick(-1:1, Nf, 'or');

rezultati3 = zeros(2*Nf+1,1);
stanje=full(RandomStateVector(2^Nf,1));

stanjet=transpose(stanje);

    for i=1:size(varijacije, 1)
        POVM = 1;
        zbir=0;

        for j=1:Nf
            t=varijacije(i,j);
            zbir=zbir+t;
            if (t==0)
                POVM = kron(POVM, E0);
            elseif (t==1)
                POVM = kron(POVM, E1);
            else
                POVM = kron(POVM, E_1);
            end
        end
        merenje=stanjet*POVM*stanje;
        rezultati3(zbir+Nf+1, 1) = rezultati3(zbir+Nf+1, 1) + merenje;
    end
    
x=(-Nf:Nf)/Nf;