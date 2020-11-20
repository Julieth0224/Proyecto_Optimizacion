function[s, si]=razon(yk, b)
% Funicion que realiza la division entre yk y b.
s = [];
si = [];
if length(b)==length(yk)
   for i=1:length(yk)
       if yk(i)>0
          m = b(i)/yk(i);
          s = [s, m];
          si = [si,i];
       end
   end
end
 
end