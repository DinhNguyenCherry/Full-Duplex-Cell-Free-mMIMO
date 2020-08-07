function [ Upsilon ] = PilotAssignment( Nm, tau, Ptr, Beta )
%PILOTASSIGNMENT Summary of this function goes here
%   Detailed explanation goes here

U = size(Beta, 2);

Beta_tilde = Nm*tau*Ptr*sum(Beta);

Upsilon = eye(tau);

if (U<=tau)
    Upsilon = Upsilon(:,1:U);
else

minHeap  = GenHeap(Beta_tilde(1:tau), Upsilon, 'min');

% nHeap = minHeap{1}
% minKey = minHeap{2}
% minParam = minHeap{3}
% minParent = minHeap{4}
% minChildren = minHeap{5}

% error('Layout is done !!!');

maxHeap = GenHeap(Beta_tilde(tau+1:U), [tau+1:U], 'max');

Upsilon = [Upsilon, zeros(size(Upsilon,1), U-tau)];

while (maxHeap{1}>0)
    Upsilon(:,maxHeap{3}(1)) = minHeap{3}(:,1);
    minHeap = RepNode(minHeap, minHeap{2}(1)+maxHeap{2}(1), minHeap{3}(:,1));
    maxHeap = RemNode(maxHeap);
%     maxHeap{1}
end
end

end

function [ ReturnHeap ] = RemNode( OrgHeap )

% Remove root node from maxHeap

CurPos = 1;
nHeap = OrgHeap{1}-1;
Key = OrgHeap{2};
Param = OrgHeap{3};
Parent = OrgHeap{4};
Children = OrgHeap{5};

while (true)
    if (sum(Children(:,CurPos)>0)==2)
        if (Key(Children(1,CurPos))>=Key(Children(2,CurPos)))
            Key(CurPos) = Key(Children(1,CurPos));
            Param(:,CurPos) = Param(:,Children(1,CurPos));
            if (sum(Children(:,Children(1,CurPos))>0)>0)
                CurPos = Children(1,CurPos);
            else
                Children(1,CurPos) = 0;
                break
            end
        else
            Key(CurPos) = Key(Children(2,CurPos));
            Param(:,CurPos) = Param(:,Children(2,CurPos));
            if (sum(Children(:,Children(2,CurPos))>0)>0)
                CurPos = Children(2,CurPos);
            else
                Children(2,CurPos) = 0;
                break
            end
        end
    elseif (sum(Children(:,CurPos)>0)==1)
        if (Children(1,CurPos)>0)
            c = 1;
        else
            c = 2;
        end
        
        Key(CurPos) = Key(Children(c,CurPos));
        Param(:,CurPos) = Param(:,Children(c,CurPos));
        if (sum(Children(:,Children(c,CurPos))>0)>0)
            CurPos = Children(c,CurPos);
        else
            Children(c,CurPos) = 0;
            break
        end
    else
        break
    end
end

ReturnHeap = {nHeap, Key, Param, Parent, Children};

end



function [ ReturnHeap ] = RepNode( OrgHeap, newKey, newParam )

% Replace root node for minHeap

CurPos = 1;
nHeap = OrgHeap{1};
Key = OrgHeap{2};
Param = OrgHeap{3};
Parent = OrgHeap{4};
Children = OrgHeap{5};

while (true)
   
    if (sum(Children(:,CurPos)>0)==2)
        if (newKey<=Key(Children(1,CurPos)) && newKey<=Key(Children(2,CurPos)))
            Key(CurPos) = newKey;
            Param(:,CurPos) = newParam;
            break
        else
            if (Key(Children(1,CurPos))<=Key(Children(2,CurPos)))
                Key(CurPos) = Key(Children(1,CurPos));
                Param(:,CurPos) = Param(:,Children(1,CurPos));
                CurPos = Children(1,CurPos);
            else
                Key(CurPos) = Key(Children(2,CurPos));
                Param(:,CurPos) = Param(:,Children(2,CurPos));
                CurPos = Children(2,CurPos);
            end
        end
    elseif (sum(Children(:,CurPos)>0)==1)
        if (Children(1,CurPos)>0)
            c = 1;
        else
            c = 2;
        end
        if (newKey<=Key(Children(c,CurPos)))
            Key(CurPos) = newKey;
            Param(:,CurPos) = newParam;
            break
        else
            Key(CurPos) = Key(Children(c,CurPos));
            Param(:,CurPos) = Param(:,Children(c,CurPos));
            CurPos = Children(c,CurPos);
        end
    else
        Key(CurPos) = newKey;
        Param(:,CurPos) = newParam;
        break
    end
        
end

ReturnHeap = {nHeap, Key, Param, Parent, Children};

end


function [ HeapTree ] = GenHeap( Keys, ParamInit, Type )

Hmin = strcmp(Type,'min');


nHeap = 0;
Key = [];
Param = [];
Parent = zeros(1, length(Keys));
Children = zeros(2, length(Keys));


for i = 1:1:length(Keys)
    Extkey = Keys(i);
    ExtPar = ParamInit(:,i);
    CurPos = 1;
    nHeap = nHeap + 1;
    while (true)
        
        if (nHeap==1)
            Key = [Key Extkey];
            Param = [Param, ExtPar];
            break
        else
            if (Hmin)
                if (Extkey<Key(CurPos))
                    tmpKey = Key(CurPos);
                    tmpPar = Param(:,CurPos);
                    Key(CurPos) = Extkey;
                    Param(:,CurPos) = ExtPar;
                    Extkey = tmpKey;
                    ExtPar = tmpPar;
                end
            else
                if (Extkey>Key(CurPos))
                    tmpKey = Key(CurPos);
                    tmpPar = Param(:,CurPos);
                    Key(CurPos) = Extkey;
                    Param(:,CurPos) = ExtPar;
                    Extkey = tmpKey;
                    ExtPar = tmpPar;
                end
            end
            
            if (Children(1,CurPos)==0)
                Parent(nHeap) = CurPos;
                Children(1,CurPos) = nHeap;
                Key(nHeap) = Extkey;
                Param(:,nHeap) = ExtPar;
                break
            elseif (Children(2,CurPos)==0)
                Parent(nHeap) = CurPos;
                Children(2,CurPos) = nHeap;
                Key(nHeap) = Extkey;
                Param(:,nHeap) = ExtPar;
                break
            else
                if (Hmin)
                    if (Key(Children(1,CurPos)) <= Key(Children(2,CurPos)) )
                        if (Extkey<=Key(Children(1,CurPos)))
                            CurPos = Children(1,CurPos);
                        else
                            CurPos = Children(2,CurPos);
                        end
                    else
                        if (Extkey<=Key(Children(2,CurPos)))
                            CurPos = Children(2,CurPos);
                        else
                            CurPos = Children(1,CurPos);
                        end
                    end
                else
                    if (Key(Children(1,CurPos)) >= Key(Children(2,CurPos)) )
                        if (Extkey>=Key(Children(1,CurPos)))
                            CurPos = Children(1,CurPos);
                        else
                            CurPos = Children(2,CurPos);
                        end
                    else
                        if (Extkey>=Key(Children(2,CurPos)))
                            CurPos = Children(2,CurPos);
                        else
                            CurPos = Children(1,CurPos);
                        end
                    end
                end
            end
            
            
        end
        
    end
    
end

HeapTree = {nHeap, Key, Param, Parent, Children};


end
