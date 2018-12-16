classdef RandomMatrix
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
           
    properties (Constant = true)
       t = 1.0e-320;
       
    end
    
    methods (Static = true) 
        
        function  V = get1dMatrix(N)
            M = eye(N);
            M(:,2:N) = M(:,2:N) + eye(N,N-1);
            M(2:N,:) = M(2:N,:) + eye(N-1,N);
            [V,D] = eig(M);
        end
        
        function interPlot(X, Y, legends, name)
           rangeX = 0.5:.001:1.3;           
           vq = interp1(X, Y, rangeX, 'pchip');
           h = figure;
           hold on;
           plot(rangeX, vq, 'LineWidth', 1.5);
           legend(legends, 2);           
           plot(X, Y, 'o', 'LineWidth', 2, 'MarkerSize',8);                    
          
           axis([0.5 1.1 0.05 1.0]);
           xlabel('\mu', 'FontSize', 28);
           ylabel('\eta', 'FontSize', 28); set(gca,'FontSize',22)
           set(gcf,'Color','w');   
           box();
           drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );
                      
           x1 = 0.75;           
           x2 = 1.03;
           
           y1 = vq(rangeX == x1,1) + 0.055;           
           y2 = vq(rangeX == x2,1) - 0.02;
                      
           ylimits1 = vq(rangeX == x1,end) - 0.03;
           ylimits2 = vq(rangeX == x2,end) + 0.03;           
           
           drawArrow([x2, x2],[y2,ylimits2],'MaxHeadSize',0.6,'Color','k','LineWidth',2); 
           drawArrow([x1, x1],[y1,ylimits1], 'MaxHeadSize',0.3, 'Color','k','LineWidth',2);
           
           hold off;           
           foutn = sprintf('%s.png',char(name));
           saveas(h, foutn, 'png');
           close(h);            
           
        end
		        
        
         function  V = get1dDisorderMatrix(N)
            M = eye(N);
            M(:,2:N) = M(:,2:N) + eye(N,N-1);
            M(2:N,:) = M(2:N,:) + eye(N-1,N);
            t = triu(rand(N)*2 - 1);
            M = M.*t + M.*t' - diag(diag(t));           
            [V,D] = eig(M);
         end
        
        function I = getIq(H,q)
           M = H;                      
           M = realpow(abs(M), 2*q); 
           Ii = sum(M);
           I = mean(Ii);
        end
        
               
        function I = getMatrixIq(N,q)
           M = RandomMatrix.get1dMatrix(N);
           I = RandomMatrix.getIq(M,q);
        end
        
        function S = getPlotVector(n0, points, q)
           N = n0;
           V = zeros(points,2);
           for i = 1:points
               V(i,1) = N;
               V(i,2) = RandomMatrix.getMatrixIq(N,q);                             
               N = N*2;               
           end
           plot(log(V(:,1)), log(V(:,2)));         
           S = V;
        end
        
        
        function V = getDisordPlotVector(n0, points, q)
           N = n0;
           V = zeros(points,2);
           for i = 1:points
               V(i,1) = N;
               M = RandomMatrix.get1dDisorderMatrix(N);               
               V(i,2) = RandomMatrix.getIq(M,q);                             
               N = N*2;               
           end
                      
        end
        
        function V = getAvgDisordPlotVector(N, points, q, times)
            V = zeros(points,2);
            V1 = zeros(points,2);
           for i = 1:times
               V1 = RandomMatrix.getDisordPlotVector(N, points, q);
               V( :, 2) = V( :,2) + V1( :,2)/times;
           end 
           V( :,1) = V1( :,1);
           plot(log(V(:,1)), log(V(:,2)))
        end
        
        function plotVectorDq(V)
            Dq = -log(V(:,2))./(log((V(:,1))));
            plot(V(:,1),Dq)
        end
        
  
        %% Figure1
        
         function I = getPreIq(H,q)
           M = H;               
           M1 = realpow(abs(M), 2*q);            
           if q == 1 
               M1 (M1 < RandomMatrix.t) = 1;
               M1 = -M1.*log(M1);           
           end               
           I = sum(sum(M1));  
           
           I = I/size(M1,2);
                        
         end
        
         function  I = getPreIvector(H,q)                  
           M1 = realpow(abs(H), 2*q);            
           if q == 1 
               M1 (M1 < RandomMatrix.t) = 1;
               M1 = -M1.*log(M1);
           end
           
           I = sum(M1);
           clear M1;
               
        end
        
        function M=getPBRM(N, b, u)
            M = triu(normrnd(0,1,N,N));			
            for i = 1:(N-1)
               for j = (i+1):N					
					M(i,j) = realpow( ( 2 * ( 1 + realpow ((sin ( pi*abs(i-j) / N ) / ( pi*b/N )), 2*u ) ) ), -0.5 ) * M(i,j);						 
					M(j,i) = M(i,j);					
               end
               M(i,i) = M(i,i)*sqrt(0.5);               
            end
            M(N,N) = M(N,N)*sqrt(0.5);
        end
        
    
       function H = addSparcity(H, a)
           N = size(H,2);
           R = triu(unifrnd(0,1,N,N));
           R = R + R' + 2*eye(N,N);                      
           H (R <= 1 - a) = 0;
           clear R;
       end
         
       function H = addSparcityMem(H, a, ident2)
           N = size(H,2);
           R = single(triu(unifrnd(0,1,N,N)));
           R = R + R' + ident2;                      
           H (R <= 1 - a) = 0;
           clear R;
       end
       
        
        %m is the 1/8 sub index of eigenstates
        function [V1] = getEigPBRM(m,M)
            N = length(M);
            [V,D] = eig(M);
			[d,ind] = sort(diag(D));
			V = V(:,ind);
            clear D d ind;
            V1 = V(:, (((N-m)/2) + 1) : ((N+m)/2));
            clear V;
        end
        
        function D = saveDATA(D, q, x)
            D(q) = x;
        end
        
          
      
      %DqVSqALL given and alfa and one u_c point 
      function [Data, eData] = getDqAlphaGraphMAT(Q,B,u,Ns,A)
            
			nMT = 18;
            nFunc = 2^(nMT - 3);             
            Data = zeros(length(Q), length(A));
            eData = zeros(length(Q), length(A));
            h = figure;
            hold on;
            for a = 1:length(A)
                sprintf('Sparsity (%f)',A(a))
                B1 = B(1); 
                Aa = A(a);
                                
                PointsQ = zeros(length(Ns), length(Q));            
                for n = 1:length(Ns)                   
                    numMatrix = 2^(nMT - Ns(n));                   
                    aux = (2^(Ns(n)-3));                   
                    N = 2^(Ns(n));                                                          
                  
                    foutA = sprintf('.\\data\\F%dB%dU%.3f_a%.2f.mat',N,B1,u,Aa);
                    
					PointsX = zeros(2, 2, 'single');					
                    if exist(foutA,'file')
                        load(foutA);  %load the data of PointsX    
                        if size(PointsX,2) ~= (nFunc)               
                            PointsX = zeros(N, nFunc, 'single');							
                        end
                    else
                        %clean corrupted data
                        PointsX = zeros(N, nFunc, 'single');     						   
                    end                     
				
                    if sum(abs(PointsX(:,nFunc))) <= 0 
                        for nM = 1: numMatrix 
                            if sum(abs(PointsX(:,nM*aux))) <= 0 
								
								ident2 = 2*eye(N,N,'single');                                
								M = single(RandomMatrix.getPBRM(N, B1, u));  
								M = RandomMatrix.addSparcityMem(M, Aa, ident2);
								%free memory
								clear ident2;
                                V0 = RandomMatrix.getEigPBRM(aux, M);
								%The central band of the RandomMatrix																
                                pIn = (nM - 1)*aux + 1;
                                pFin = (nM)*aux;
                                PointsX(:,pIn:pFin) = V0;
								clear V0;
								clear M;                            
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						save(foutA, 'PointsX','-v7.3');						
                    end
                     
                    for q = 1:length(Q)                                        
						
                        V = RandomMatrix.getPreIvector(PointsX, Q(q)); 						
                        if Q(q) ~= 1						
						 V = log(V);
                        end
                        PointsQ(n,q) = mean(V);                        
                    end                    
                    line = sprintf('(%d) N %d nMatrix %d  Eig %d  X %d',n, N, numMatrix, aux, size(PointsX,2));
					disp(line);
					clear PointsX;					
                end            
              
                for q = 1:length(Q)
				
                    Dq = zeros(1,2);
                    S.normr = 10;
                    if Q(q) == 1
                        [Dq, S] = polyfit(log(realpow(2,Ns)), PointsQ(:,q)',1);                         
                    else
                        [Dq, S] = polyfit(log(realpow(2,Ns)), PointsQ(:,q)'./-(Q(q)-1),1); 
                   
                    end
                    Data(q,a) = Dq(1);
                    eData(q,a) = S.normr;
                   
                end
                
				lineName = sprintf('\\alpha = %.1f', A(a));
                errorbar(Q,Data(:,a),eData(:,a), 'DisplayName',char(lineName), 'LineWidth', 2);				
				disp('     DATA       eDATA');
				disp([Data(:,a), eData(:,a)]);
				
            end 
           xlabel('q', 'FontSize', 24);
           ylabel('D_q', 'FontSize', 24); set(gca,'FontSize',18)
           set(gcf,'Color','w');   
           box on;
		   legend('show');		   
           hold off;                  
           		   
		   % saving figure
		   saveas(h, 'DqVSq.png', 'png');		   
           close(h);
            
        end  
      
         
        %only for one u 
        function [Data, eData] = getDqGraph(Q,B,u,Ns)
            
            nMT = 18;
            Data = zeros(length(B), length(Q));
            eData = zeros(length(B), length(Q));
            
            Points = zeros(length(Ns), length(Q));
            
            for b = 1:length(B)
                sprintf('BAND (%f)',B(b))
                Points = zeros(length(Ns), length(Q));
                for n = 1:length(Ns)                   
                    numMatrix = 2^(nMT - Ns(n));                   
                    aux = ( 2^(Ns(n)-2));
                    %X = (numMatrix * aux);
                    N = 2^(Ns(n));
                    for nM = 1: numMatrix                       
                        M = RandomMatrix.getPBRM(N, B(b), u);                       
                        V = RandomMatrix.getEigPBRM(aux, M);                        
                        for q = 1:length(Q)
                            Points(n,q) = Points(n,q) + RandomMatrix.getPreIq(V, Q(q))/numMatrix;                    
                        end
                    end                  
                                          
                    line = sprintf('(%d) N %d nMatrix %d  Eig %d  X %d',n, N, numMatrix, aux, numMatrix*size(V,2))                      
                end            
              
                for q = 1:length(Q)
                    Dq = zeros(1,2);
                    S.normr = 10;
                    if Q(q) == 1
                        [Dq, S] = polyfit(log(realpow(2,Ns)), Points(:,q)',1);                         
                    else                       
                        [Dq, S] = polyfit(log(realpow(2,Ns)), log(Points(:,q))'./-(Q(q)-1),1);
                        %[Dq, S] = polyfit(log(realpow(2,Ns)), Points(:,q)'./-(Q(q)-1),1);
                    end
                    %Data(b,q) = Dq(1);
                    Data(b,:) = RandomMatrix.saveDATA(Data(b,:), q, Dq(1));                   
                    eData(b,:) = RandomMatrix.saveDATA(eData(b,:), q, S.normr);                                       
                end
                
            end           
        end
        
        function [col, PointsL, PointsX] = getColPoints(PointsL, PointsX, N)
            col = find(PointsL == N);
            if isempty(col) == 1
                col = size(PointsL,2) + 1;                
				PointsL = [PointsL, N];
				PointsX = [PointsX, zeros(size(PointsX,1), 1, 'single')];
            end
            return;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [PointsVar] = getEtaUGraph(Q,B,U,Ns,a)
    
            NT = maxNumCompThreads;
            maxNumCompThreads(NT);
            nMT = 18;
            
            name = sprintf('sB%da%.2f.txt', B(1), a);
            PointsVar = zeros(length(U), length(Ns), 'single');  
            nFunc = 2^(nMT - 3); 
                        
            for u = 1:length(U)
                                
                lineStat = sprintf('U (%f)',U(u));
                disp(lineStat);                
                foutF = sprintf('.\\data\\sB%da%.2fU%.3f.mat', B(1), a, U(u));
                
                if exist(foutF,'file')
					load(foutF);
					if size(PointsX,1) ~= (nFunc)
						PointsX = zeros(nFunc, length(Ns), 'single');
						PointsL = single(2.^(Ns));
					end
				else
					PointsX = zeros(nFunc, length(Ns), 'single');
					PointsL = single(2.^(Ns));
                end
                
                for n = 1:length(Ns)
                                        
                    numMatrix = 2^(nMT - Ns(n));                   
                    aux = ( 2^(Ns(n)-3));
                    N = 2^(Ns(n));                  
                    V = [];                    
                    Q1 = Q(1);
                    B1 = B(1);
                    Uu = U(u);
                    [n1,PointsL,PointsX] = RandomMatrix.getColPoints(PointsL, PointsX, N);					
                    
                    if  sum(PointsX(:,n1)) <= 0 
                        parfor nM = 1: numMatrix  
							ident2 = 2*eye(N,N,'single'); 
                            M = single(RandomMatrix.getPBRM(N, B1, Uu));
                            M = RandomMatrix.addSparcityMem(M, a, ident2);							
                            V0 = RandomMatrix.getEigPBRM(aux, M);                           
                            V = [V, RandomMatrix.getPreIvector(V0, Q1)];                        
                        end
                        PointsX(:,n1) = V';
						save(foutF, 'PointsX','PointsL','-v7.3');
                        
                    end
                    
                    Var1 = var(PointsX(:,n1),1);
                    deltaXi = sqrt(Var1);
                    meanXi = mean(PointsX(:,n1));                    
                    lineStat = sprintf('Size %d  numMatrix %d \n [ var=%f  mean=%f X = %f ]', N,  numMatrix, Var1, meanXi, deltaXi/meanXi);
                    disp(lineStat);                                        
                    PointsVar(u,:) = RandomMatrix.saveDATA(PointsVar(u,:), n, deltaXi/meanXi);                 
                    nPoints = [(U)', PointsVar];                    
                    dlmwrite(name, nPoints);                        
                end               
            end
                        
            legends =  char(sprintf('%d',2^(Ns(1))));
            for n = 2:length(Ns)
               legends = char(legends, sprintf('%d',2^(Ns(n)) )); 
            end            
            RandomMatrix.interPlot(nPoints(:,1), nPoints(:,2:size(nPoints,2)), legends, 'sBaFigure');            
        end

        
      
            
    end
    
end

