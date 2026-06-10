%Diego Oliva, Erik Cuevas, Gonzalo Pajares, Daniel Zaldivar y Marco Perez-Cisneros
%Multilevel Thresholding Segmentation Based on Harmony Search Optimization
%Universidad Complutense de Madrid / Universidad de Guadalajara

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The algorithm was published in:
%Diego Oliva, Erik Cuevas, Gonzalo Pajares, Daniel Zaldivar, and Marco Perez-Cisneros, 
%“Multilevel Thresholding Segmentation Based on Harmony Search Optimization,” 
%Journal of Applied Mathematics, vol. 2013, 
%Article ID 575414, 24 pages, 2013. doi:10.1155/2013/575414
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fitR, fitBestR, fitG, fitBestG, fitB, fitBestB] = Otsu(I, N, Lmax, level, xR, probR)
%Metodo de Otsu

%Evalua poblaciones xR, xG, xB, en la funcion objetivo para obtener
%fitR, fitG, fitB, dependiendo si la imagen es RGB o escala de grises
for j = 1:N

        %RGB image
        fitR(j) = sum(probR(1:xR(j,1))) * (sum((1:xR(j,1)) .* probR(1:xR(j,1)) / sum(probR(1:xR(j,1)))) - sum((1:Lmax) .* probR(1:Lmax))) ^ 2;
        for jlevel = 2:level - 1
            fitR(j) = fitR(j) + sum(probR(xR(j,jlevel - 1) + 1:xR(j,jlevel))) * (sum((xR(j,jlevel - 1) + 1:xR(j,jlevel)) .* probR(xR(j,jlevel - 1) + 1:xR(j,jlevel)) / sum(probR(xR(j,jlevel - 1) + 1:xR(j,jlevel)))) - sum((1:Lmax) .* probR(1:Lmax))) ^ 2;
        end
        fitR(j) = fitR(j) + sum(probR(xR(j,level-1) + 1:Lmax)) * (sum((xR(j,level - 1) + 1:Lmax) .* probR(xR(j,level - 1) + 1:Lmax) / sum(probR(xR(j,level - 1) + 1:Lmax))) - sum((1:Lmax) .* probR(1:Lmax))) ^ 2;

        fitBestR(j) = fitR(j);
        
%         fitG(j) = sum(probG(1:xG(j,1))) * (sum((1:xG(j,1)) .* probG(1:xG(j,1)) / sum(probG(1:xG(j,1)))) - sum((1:Lmax) .* probG(1:Lmax))) ^ 2;
%         for jlevel = 2:level - 1
%             fitG(j) = fitG(j) + sum(probG(xG(j,jlevel - 1) + 1:xG(j,jlevel))) * (sum((xG(j,jlevel - 1) + 1:xG(j,jlevel)) .* probG(xG(j,jlevel - 1) + 1:xG(j,jlevel)) / sum(probG(xG(j,jlevel - 1) + 1:xG(j,jlevel)))) - sum((1:Lmax) .* probG(1:Lmax))) ^ 2;
%         end
%         fitG(j) = fitG(j) + sum(probG(xG(j,level - 1) + 1:Lmax)) * (sum((xG(j,level-1) + 1:Lmax) .* probG(xG(j,level - 1) + 1:Lmax) / sum(probG(xG(j,level - 1) + 1:Lmax))) - sum((1:Lmax) .* probG(1:Lmax))) ^ 2;
%         fitBestG(j) = fitG(j);
        
%         fitB(j) = sum(probB(1:xB(j,1))) * (sum((1:xB(j,1)) .* probB(1:xB(j,1)) / sum(probB(1:xB(j,1)))) - sum((1:Lmax) .* probB(1:Lmax))) ^ 2;
%         for jlevel = 2:level - 1
%             fitB(j) = fitB(j) + sum(probB(xB(j,jlevel - 1) + 1:xB(j,jlevel))) * (sum((xB(j,jlevel - 1) + 1:xB(j,jlevel)) .* probB(xB(j,jlevel - 1) + 1:xB(j,jlevel)) / sum(probB(xB(j,jlevel - 1) + 1:xB(j,jlevel)))) - sum((1:Lmax) .* probB(1:Lmax))) ^ 2;
%         end
%         fitB(j) = fitB(j) + sum(probB(xB(j,level - 1) + 1:Lmax)) * (sum((xB(j,level - 1) + 1:Lmax) .* probB(xB(j,level - 1) + 1:Lmax) / sum(probB(xB(j,level - 1) + 1:Lmax))) - sum((1:Lmax) .* probB(1:Lmax))) ^ 2;
%         fitBestB(j) = fitB(j);
    
end
  
      %Imagen escala de Grises
        fitR = fitR';
        fitBestR = fitBestR';
