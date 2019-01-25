
function [gbest,gbestval,C,fitcount]= mqhao_with_subgroup2(Dimension,Particle_Number,g,gr,Max_Gen,VRmin,VRmax,varargin)
%对种群进行分组进化比对

DIM=Dimension;
minDomain=VRmin;
maxDomain=VRmax;
maxFE=Max_Gen;

repeat=1;
C=zeros(1,310000);


group=g;
gr=gr;%每组数目
% groupNum=gr*group;
optimalNum=gr*group;%k值的大小
%取前alfa中的最优解作为下一代的向标
alfa=0.5;

%------参数定义区域结束------
%------多次重复计算所需参数------
gbestV=zeros(1,repeat);  gfe=zeros(1,repeat); tot_time1 = zeros(1,repeat);
%------多次重复计算------
for kk=1:repeat
    tic;
    %------参数初始化区域开始------
    funcV=zeros(1,optimalNum);  %定义保存k个采样点的函数值，初始化为0
    samplePos=zeros(optimalNum,DIM);%采样点矩阵
    sigma=maxDomain-minDomain; %初始尺度为目标函数定义域的大小
    stdPre=zeros(1,DIM);%存储每一个维度上的标准差
    stdNow=zeros(1,DIM);
    %------计算k个采样位置的目标函数值------
    optimalSolution=unifrnd(minDomain,maxDomain,optimalNum,DIM);%定义k个DIM维采样点的坐标，并初始化
    stdPre=std(optimalSolution,1,1) ;%计算初始时k个采样点的标准差，按行求
    w=1;% function evolution times
    
    VRmin=repmat(minDomain,gr,DIM);
    VRmax=repmat(maxDomain,gr,DIM);
    
   
    for k=1:optimalNum  %求最优解函数值
        funcV(k)=func(optimalSolution(k,:),DIM,varargin{:});
        
    end
     
%     optimalSolution=optimalSolution(index_sort,:);
    %更新prefix_optimalSolution 和计算prefix_optimalSolution mean
    [~,index_sort]=sort(funcV);
    prefix_optimalSolution =optimalSolution(index_sort(1:ceil(alfa*optimalNum)),:);
    first_prefix_optimalSolution_mean=1/(ceil(alfa*optimalNum)).*(sum(prefix_optimalSolution));
    
    mea=first_prefix_optimalSolution_mean;
    C(w)=min(funcV);
    csigma = ones(1,DIM)*(sigma);
    covv = diag(csigma.^2);
%     covv=cov(prefix_optimalSolution);
    while 1   % M iteration begin
        if w>maxFE
            break;
        end
        while 1  %QHO iteration begin
            if w>maxFE
                break;
            end
            while 1  %稳定性收敛迭代开始
                
                if w>maxFE
                    break;
                end
                %                change_flag=0; % op_solution更新判断标志
                %                  for k=1:optimalNum
                for g=1:group
                    while 1
                        premean=mea;
                         %计算prefix_optimalSolution协方差，使用上一次的均值计算,算出来的协方差矩阵不正定
                         %使用mvnrnd会报错
                        covv=coov(prefix_optimalSolution,premean);
                        change_flag=0; % op_solution更新判断标志
                      
                        samplePos = mvnrnd(mean(prefix_optimalSolution,1),covv,gr);
                        samplePos=(samplePos>VRmax).*VRmax+(samplePos<=VRmax).*samplePos;
                        samplePos=(samplePos<VRmin).*VRmin+(samplePos>=VRmin).*samplePos;
                        
                        %                     for  i=1:optimalNum
                        for i=1:gr%1:groupNum
                            sampleValue=func(samplePos(i,:),DIM,varargin{:});%求第i个采样点的函数值
                            w=w+1;
                            C(w)=min(funcV);
                            if sampleValue<funcV((((g-1)*gr+i)))     %如果采样点值小于当前点函数值，则替换,
                                funcV((((g-1)*gr+i)))=sampleValue;
                                optimalSolution((((g-1)*gr+i)),:)=samplePos(i,:);
                                change_flag=1;
                            end
                        end
                        %                  end
                        if(change_flag==0)
                            
                            break;
                        end
                    end %稳定性收敛迭代结束
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
            [v_max,index_max]=max(funcV);%取得最大值的序号index_max
            optimalSolution(index_max,:)=nom;%用平均坐标替换最大值对应坐标
            funcV(index_max)=func(nom,DIM,varargin{:});%
            w=w+1;
            C(w)=min(funcV);
            stdPre=std(optimalSolution,1,1);  %新解标准差

            %更新prefix_optimalSolution
             [~,index1]=sort(funcV);
            prefix_optimalSolution =optimalSolution( index1(1:ceil(alfa*optimalNum)),:);
            %计算本代均值，供下次计算方差使用
            mea=mean(optimalSolution );

            %本代自己纯协方差
%          covv=cov(prefix_optimalSolution);


            %------均值替换
            if w>maxFE
                break;
            end
            % ------能级下降结束------
            if max(stdPre)<sigma
                break;
            end
        end % QHO iteration end
        sigma=sigma/2.0;
        

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

