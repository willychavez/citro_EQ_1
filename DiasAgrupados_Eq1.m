clear all
clc
raiz = pwd;
cd('Link to Equipamento_1')
ym= dlmread("YMedios_Eq1.txt");
cd('Agrupados_Eq1/09_22_2017');

%% CARREGA DADOS

files=dir;

for k=1:(length(files)-2)
cd(files(k+2).name);
arquivos=dir('*.txt'); % identifica os arquivos .txt na pasta
  for j=1:length(arquivos)
    data=importdata(arquivos(j).name,'\t',1);     
    YYx{k}(:,j)=data.data(:,2); % Salva os espectros na matriz 
  end   
cd ..
end
xx=data.data(:,1); % Salva os comprimentos de onda no vetor

clear k i j

%Corte
 ind1 = find (xx > 470 & xx < 770);
for i=1:length(YYx)
  for j=1:size(YYx{i},2)
    YYr{i}(:,j) = YYx{i}(ind1,j);
  end
end 
xxv=xx(ind1,:);
%figure; plot(xxv, YYr{3});  

clear k i j

%filtro boxcar=11
fboxcar=10;
for k=1:length(YYx)
  for j=1:size(YYr{k},2)
    for i=(1+fboxcar/2):(size(YYr{k},1)-fboxcar/2)
     YYv{k}(i,j) = mean(YYr{k}((i-fboxcar/2):(i+fboxcar/2),j));
    end
  end
end 

clear k i j  

%% CORTA O ESPECTRO

%Equipamento 1
ind = find (xxv > 494.3046 & xxv < 758.7687);
for i=1:length(YYv)
    for j=1:size(YYv{i},2)
        YY{i}(:,j) = YYv{i}(ind,j);
    end
end
x=xxv(ind,:);

clear j i

%Retira o offset
for i=1:length(YY)
  for j=1:size(YY{i},2)
    m=min(YY{i}(:,j));
    YYoff{i}(:,j)=YY{i}(:,j)-m;
  end
end

clear j i

%Produto Escalar
A1 = [];
A2 = [];
A3 = [];

for j=1:size(ym,2)
norma_ym = norm(ym(:,j));
  for i=1:size(YYoff{j},2)
    norma_y = norm(YYoff{j}(:,i));
    costheta = dot(YYoff{j}(:,i),ym(:,j))/(norma_y * norma_ym);
    if j==1 &&  costheta < 0.992
      A1 = [A1;i];
    endif
    if j==2 &&  costheta < 0.995
      A2 = [A2;i];
    endif
    if j==3 &&  costheta < 0.989
      A3 = [A3;i];
    endif
  end
end

pause
clear k j

%identifica os espectros das folhas
cd([raiz '/Link to Equipamento_1/Numeros_folhas_W']);
files1=dir;
b1=[];
b2=[];
for k=1:(length(files)-2)
  cd(files1(k+2).name);
  arquivos1=dir('*.txt'); 
  for j=1:1%length(arquivos1)
    m=importdata(arquivos1(j).name,'\t');
    b1=[b1;m];
  end
  n_folhas{k}(:,1)=b1;
  b1=[];
  cd ..
end

clear k j

for k=1:(length(files1)-2)
  for i=1:length(n_folhas{k}())
     bi = [i*ones(n_folhas{k}(i,1),1), [1:n_folhas{k}(i,1)]'];
     b2=[b2;bi];
  end
  folhas{k}(:,1:2)=b2;
  b2=[];
end

%Remove outliers
clear s j a

for s=1:length(A1)
  j = s - 1;
  a = A1(s) - j;
  YYoff{1}(:,a)=[];
  folhas{1}(a,:)=[];
end  

clear s j a

for s=1:length(A2)
  j = s - 1;
  a = A2(s) - j;
  YYoff{2}(:,a)=[];
  folhas{2}(a,:)=[];
end  


clear s j a 

for s=1:length(A3)
  j = s - 1;
  a = A3(s) - j;
  YYoff{3}(:,a)=[];
  folhas{3}(a,:)=[];
end


clear j i

%% NORMALIZA PELA AREA
for i=1:length(YYoff)
   for j=1:size(YYoff{i},2)
      %Calcula area embaixo da reta
      int{i}=trapz(x,YYoff{i});
      % Divido cada valor de intensidade do pico pela area abaixo da RETA
      YYcorrigido{i}(:,j)=YYoff{i}(:,j)/int{i}(1,j);
   end
end

cd([raiz '/Link to Equipamento_1/Medidas']);

Espectros_1=[x,YYcorrigido{1}];
Espectros_2=[x,YYcorrigido{2}];
Espectros_3=[x,YYcorrigido{3}];
dlmwrite ('ASSIN_09_22_2017_Eq1.dat', Espectros_1, 'delimiter','\t', 'precision',4)
dlmwrite ('SADIA_09_22_2017_Eq1.dat', Espectros_2, 'delimiter','\t', 'precision',4)
dlmwrite ('SINT_09_22_2017_Eq1.dat', Espectros_3, 'delimiter','\t', 'precision',4)

figure;plot(x,YYcorrigido{1})
figure;plot(x,YYcorrigido{2})
figure;plot(x,YYcorrigido{3})

cd([raiz '/Link to Equipamento_1/weka']);

%% ARFF
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%nome com que o arquivo arff deve ser salvo!
%cd ..

fid = (fopen('09_22_2017_Eq1.arff','w'));
fprintf(fid,'@RELATION Diagnostico \n');
fprintf(fid,'\n');
%o primeiro atributo sera o nome da folha
fprintf(fid,'@ATTRIBUTE PLANTA string \n');
fprintf(fid,'@ATTRIBUTE FOLHA string \n');
fprintf(fid,'\n');
for i=1:size(x,1)
    fprintf(fid,'@ATTRIBUTE %s NUMERIC \n',num2str(x(i)));
end

fprintf(fid,'@ATTRIBUTE class {Assintomatica,Sadia,Sintomatica} \n');

fprintf(fid,'\n');
fprintf(fid,'@DATA \n');

    
%fprintf(fid,'%d,',k);    
for l=1:length(YYcorrigido{1}(1,:))
  fprintf(fid,'%i,',folhas{1}(l,1:end));
  for j=1:size(x,1)
    fprintf(fid,'%0.5f,',YYcorrigido{1}(j,l));
  end
  fprintf(fid,'Assintomatica');
  fprintf(fid,'\n');
end                            

%fprintf(fid,'%d,',k);    
for l=1:length(YYcorrigido{2}(1,:))
  fprintf(fid,'%i,',folhas{2}(l,1:end));
  for j=1:size(x,1)
    fprintf(fid,'%0.5f,',YYcorrigido{2}(j,l));
  end
  fprintf(fid,'Sadia');
  fprintf(fid,'\n');
end

%fprintf(fid,'%d,',k);    
for l=1:length(YYcorrigido{3}(1,:))
  fprintf(fid,'%i,',folhas{3}(l,1:end));
  for j=1:size(x,1)
    fprintf(fid,'%0.5f,',YYcorrigido{3}(j,l));
  end
  fprintf(fid,'Sintomatica');
  fprintf(fid,'\n');
end 
fclose(fid);    

cd(raiz)

