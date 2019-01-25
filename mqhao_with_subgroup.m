
function [gbest,gbestval,C,fitcount]= mqhao_with_subgroup(Dimension,Particle_Number,g,gr,Max_Gen,VRmin,VRmax,varargin)
%����Ⱥ���з�������ȶ�

DIM=Dimension;
minDomain=VRmin;
maxDomain=VRmax;
maxFE=Max_Gen;

repeat=1;
C=zeros(1,310000);


group=g;
gr=gr;%ÿ����Ŀ
% groupNum=gr*group;
optimalNum=gr*group;%kֵ�Ĵ�С

%------���������������------
%------����ظ������������------
gbestV=zeros(1,repeat);  gfe=zeros(1,repeat); tot_time1 = zeros(1,repeat);
%------����ظ�����------
for kk=1:repeat
    tic;
    %------������ʼ������ʼ------
    funcV=zeros(1,optimalNum);  %���屣��k��������ĺ���ֵ����ʼ��Ϊ0
    samplePos=zeros(optimalNum,DIM);%���������
    sigma=maxDomain-minDomain; %��ʼ�߶�ΪĿ�꺯��������Ĵ�С
    stdPre=zeros(1,DIM);%�洢ÿһ��ά���ϵı�׼��
    stdNow=zeros(1,DIM);
    %------����k������λ�õ�Ŀ�꺯��ֵ------
    optimalSolution=unifrnd(minDomain,maxDomain,optimalNum,DIM);%����k��DIMά����������꣬����ʼ��
    stdPre=std(optimalSolution,1,1) ;%�����ʼʱk��������ı�׼�������
    w=1;% function evolution times
    
    VRmin=repmat(minDomain,gr,DIM);
    VRmax=repmat(maxDomain,gr,DIM);
    for k=1:optimalNum  %�����Ž⺯��ֵ
        funcV(k)=func(optimalSolution(k,:),DIM,varargin{:});
        
    end
    C(w)=min(funcV);
    csigma = ones(1,DIM)*(sigma);
    covv = diag(csigma.^2);
%     covv=cov(optimalSolution);
    while 1   % M iteration begin
        if w>maxFE
            break;
        end
        while 1  %QHO iteration begin
            if w>maxFE
                break;
            end
            while 1  %�ȶ�������������ʼ
                
                if w>maxFE
                    break;
                end
                %                change_flag=0; % op_solution�����жϱ�־
                %                  for k=1:optimalNum
                for g=1:group
                    while 1
                        change_flag=0; % op_solution�����жϱ�־
                        samplePos = mvnrnd(mean(optimalSolution((((g-1)*gr+1)):(g*gr),:),1),covv,gr);
                        
                        samplePos=(samplePos>VRmax).*VRmax+(samplePos<=VRmax).*samplePos;
                        samplePos=(samplePos<VRmin).*VRmin+(samplePos>=VRmin).*samplePos;
                        
                        %                     for  i=1:optimalNum
                        for i=1:gr%1:groupNum
                            sampleValue=func(samplePos(i,:),DIM,varargin{:});%���i��������ĺ���ֵ
                            w=w+1;
                            C(w)=min(funcV);
                            if sampleValue<funcV((((g-1)*gr+i)))     %���������ֵС�ڵ�ǰ�㺯��ֵ�����滻,
                                funcV((((g-1)*gr+i)))=sampleValue;
                                optimalSolution((((g-1)*gr+i)),:)=samplePos(i,:);
                                change_flag=1;
                            end
                        end
                        %                  end
                        if(change_flag==0)
                            
                            break;
                        end
                    end %�ȶ���������������
                end
                break;
            end
            
            
            swarmNum=optimalNum;
            a=1:swarmNum;
            alpha=sum(log((swarmNum+1)./(2*a)));
            weights=(log((swarmNum+1)./(2.*a)))/alpha;
            [v,vIndex] = sort(funcV);
            nm=zeros(swarmNum,DIM);
            for ii=1:swarmNum
                nm(ii,:)=nm(ii,:)+weights(ii)'.*optimalSolution(vIndex(ii),:);
            end
            nom=sum(nm,1);
            [v_max,index_max]=max(funcV);%ȡ�����ֵ�����index_max
            optimalSolution(index_max,:)=nom;%��ƽ�������滻���ֵ��Ӧ����
            funcV(index_max)=func(nom,DIM,varargin{:});%
            w=w+1;
            C(w)=min(funcV);
            stdPre=std(optimalSolution,1,1);  %�½��׼��
            
            %------��ֵ�滻
            if w>maxFE
                break;
            end
            % ------�ܼ��½�����------
            if max(stdPre)<sigma
                break;
            end
        end % QHO iteration end
        sigma=sigma/2.0;
        
        csigma = ones(1,DIM)*(sigma);
        covv = diag(csigma.^2);
%         covv=cov(optimalSolution);
        %        if sigma<=ACCURACY
        %            break;
        %        end
        if w>maxFE
            break;
        end
    end % M iteration end
    %optimalSolution
    tot_time1(kk) = toc;
    
    [global_min,index]=min(funcV);
    gbest=min(funcV);
    gbestval=min(funcV);
    fitcount=w;
end % MQHOA end

