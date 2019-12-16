% function opt = bisection_grad()
%     grad = diff(p, u);
%     a = 0;
%     b = 0.1;
%     l = 0.001;
%     n = ceil(-log2( l/(b-a)));
% 
%     for k = 1:n
%         x0 = (a + b)/2;
%         if grad(x0) == 0
%             a=x0;
%             b=x0;
%             break;
%         elseif grad(x0) > 0
%             a = a;
%             b = x0;
%         else
%             a = x0;
%             b = b;
%         end
%     end
%     opt = x0;
% end
