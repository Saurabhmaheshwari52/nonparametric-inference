function plot_density(a1,a2,a3,ind, kind)
    p = kind;
    val = normpdf(sort(a1(:,ind)), mean(a1(:,ind)), std(a1(:,ind)));
    plot(sort(a1(:,ind)), val,'-or')
    grid on
    hold on
    val = normpdf(sort(a2(:,ind)), mean(a2(:,ind)), std(a2(:,ind)));
    plot(sort(a2(:,ind)), val,'-ob')
    grid on
    hold on
    val = normpdf(sort(a3(:,ind)), mean(a3(:,ind)), std(a3(:,ind)));
    plot(sort(a3(:,ind)), val,'-og')
    grid on
    if strcmp(p,'normal')
        loc = 1;
%         legend('scale = 0.2','scale = 0.3','scale = 0.4')
    else if strcmp(p,'lognormal')
         loc = -0.1;
%          legend('scale = 0.2','scale = 0.3','scale = 0.4')
    else if strcmp(p,'weibull')
         loc = 0.5;
%          legend('scale = 1.2','scale = 1.35','scale = 1.5')
        end
        end
    end
    if ind > 3
            title(strcat('x', num2str(ind-3),'scale ', '(True loc = ', num2str(loc),')'))
            
    else
            title(strcat('x', num2str(ind),'loc ', '(True loc = ', num2str(loc),')'))
    end    
    
end

