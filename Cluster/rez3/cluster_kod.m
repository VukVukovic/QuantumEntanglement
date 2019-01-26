addpath(genpath('/home/polaznik17/upload'));

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

Nf=13;
    varijacije=pick(-1:1, Nf, 'or'); 
            qubit1=1;
            for j=1:Nf
                qubit1=Tensor(qubit1, abs(RandomStateVector(2,1)));
            end
            qubit2=full(GHZState(2,Nf));
            qubit3=full(DickeState(Nf,2));
            qubit4=full(DickeState(Nf,3));
            qubit5=full(WState(Nf));
            qubit6=full(RandomStateVector(2^Nf,1));
            
        qubitt1=transpose(qubit1);
        qubitt2=transpose(qubit2);
        qubitt3=transpose(qubit3);
        qubitt4=transpose(qubit4);
        qubitt5=transpose(qubit5);
        qubitt6=transpose(qubit6);

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
        rezultati1(i, 1) = qubitt1*POVM*qubit1;
        rezultati2(i, 1) = qubitt2*POVM*qubit2;
        rezultati3(i, 1) = qubitt3*POVM*qubit3;
        rezultati4(i, 1) = qubitt4*POVM*qubit4;
        rezultati5(i, 1) = qubitt5*POVM*qubit5;
        rezultati6(i, 1) = qubitt6*POVM*qubit6;
        zbirovi(i, 1) = zbir;
        end

        rezultatifinal = zeros(2*Nf+1,1);
        for i=1:size(varijacije, 1)
            rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati1(i);
        end
        
        filename=strcat('rezultati1.mat');
        save(filename, 'rezultatifinal');
        
        rezultatifinal = zeros(2*Nf+1,1);
        for i=1:size(varijacije, 1)
            rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati2(i);
        end
        
        filename=strcat('rezultati2.mat');
        save(filename, 'rezultatifinal');
        
        rezultatifinal = zeros(2*Nf+1,1);
        for i=1:size(varijacije, 1)
            rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati3(i);
        end
        
        filename=strcat('rezultati3.mat');
        save(filename, 'rezultatifinal');
        
        rezultatifinal = zeros(2*Nf+1,1);
        for i=1:size(varijacije, 1)
            rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati4(i);
        end
        
        filename=strcat('rezultati4.mat');
        save(filename, 'rezultatifinal');
        
        rezultatifinal = zeros(2*Nf+1,1);
        for i=1:size(varijacije, 1)
            rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati5(i);
        end
        
        filename=strcat('rezultati5.mat');
        save(filename, 'rezultatifinal');
        
        rezultatifinal = zeros(2*Nf+1,1);
        for i=1:size(varijacije, 1)
            rezultatifinal(zbirovi(i,1)+Nf+1)=rezultatifinal(zbirovi(i,1)+Nf+1)+rezultati6(i);
        end
        
        filename=strcat('rezultati6.mat');
        save(filename, 'rezultatifinal');