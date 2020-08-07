%% Small-scale Fading & Channel Model

%% DL and UL

H_A2D = cell(K, M);
H_A2U = cell(M, L);

H_DL = [];
H_UL = [];

for m = 1:1:M
    
    temp_DL = [];
    for k = 1:1:K
        H_A2D{k,m} = CreateSmallScaleFading( 1, Nm );
        temp_DL = [temp_DL, H_A2D{k,m}'*D_A2D(m,k)];
    end
    H_DL = [H_DL; temp_DL];
    
    temp_UL = [];
    for l = 1:1:L
        H_A2U{m,l} = CreateSmallScaleFading( Nm, 1 );
        temp_UL = [temp_UL, H_A2U{m,l}*D_A2U(m,l)];
    end
    H_UL = [H_UL; temp_UL];
    
end

H_DL = H_DL';


%% A2A and SI channels

G_AA = cell(M, M);
G_AA_All = [];

for m = 1:1:M
    temp_AA = [];
    for mp = 1:1:M
        if (m==mp)
            v = sqrt(10); s = 1; 
            dist = {'Rice', v, s}; % K_factor = v^2/(2*s^2) = 5
            G_AA{m, mp} = CreateSmallScaleFading( Nm, Nm, dist );
%             G_AA{m, mp} = CreateSmallScaleFading( Nm, Nm );
        else
            G_AA{m, mp} = CreateSmallScaleFading( Nm, Nm ) * D_A2A(mp, m);
        end
        temp_AA = [temp_AA; G_AA{m, mp}];
    end
    G_AA_All = [G_AA_All, temp_AA];
end


G_CCI = CreateSmallScaleFading( K, L ) .* D_U2D(1:L,1:K)';



