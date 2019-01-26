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

Nf=8;
J=6;

varijacije = pick(-1:1, Nf, 'or');

qubit=full(DickeState(Nf, J));
qubitt=transpose(qubit);

parfor i=1:size(varijacije, 1)
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
    merenje=qubitt*POVM*qubit;
    rezultati(i, 1) = merenje;
    zbirovi(i, 1) = zbir;
end

rezultatifinal = zeros(2*Nf+1,1);
for i=1:size(varijacije, 1)
    rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati(i);
end

pr = sum(rezultatifinal);
X=(-Nf:Nf)/Nf;
plot(X, rezultatifinal);