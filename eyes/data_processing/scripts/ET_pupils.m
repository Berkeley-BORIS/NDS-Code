clear all; close all;
figure(); hold on;

file1 = {'../data/gesdrdp1/gesdrdp1_fixation_points_raw_pupil.mat', ...
    '../data/bwsdrdp6/bwsdrdp6_fixation_points_raw_pupil.mat', ...
    '../data/bwsdrdp8/bwsdrdp8_fixation_points_raw_pupil.mat'};
file2 = {'/Users/natdispstats/Documents/eyes/data_processing/data/gesdrdp1/gesdrdp1_proc.mat', ...
    '/Users/natdispstats/Documents/eyes/data_processing/data/bwsdrdp6/bwsdrdp6_proc.mat', ...
    '/Users/natdispstats/Documents/eyes/data_processing/data/bwsdrdp8/bwsdrdp8_proc.mat'};
for f = 1:3
    
    load(file1{f});
    load(file2{f});
    %load(['../data/gesdrdp1/gesdrdp1_fixation_points_raw_pupil.mat'])
    %load(['/Users/natdispstats/Documents/eyes/data_processing/data/gesdrdp1/gesdrdp1_proc.mat'])
    
    all50 = [radials_50_1 ; radials_50_2 ; radials_50_3 ; radials_50_4];
    all100 = [radials_100_1 ; radials_100_2 ; radials_100_3 ; radials_100_4];
    all450 = [radials_450_1 ; radials_450_2 ; radials_450_3 ; radials_450_4];
    alltargets = [all50 ; all100 ; radials_200_3 ; all450];
    allwalks = [walk_1 ; walk_2 ; walk_3];
    
    % figure(); hold on;
    %
    % display (['50cm means']);
    % display (['rep 1    Left: ' num2str(mean(radials_50_1(:,17))) '   Right: ' num2str(mean(radials_50_1(:,18)))]);
    % display (['rep 2    Left: ' num2str(mean(radials_50_2(:,17))) '   Right: ' num2str(mean(radials_50_2(:,18)))]);
    % display (['rep 3    Left: ' num2str(mean(radials_50_3(:,17))) '   Right: ' num2str(mean(radials_50_3(:,18)))]);
    % display (['rep 4    Left: ' num2str(mean(radials_50_4(:,17))) '   Right: ' num2str(mean(radials_50_4(:,18)))]);
    % display (['ALL      Left: ' num2str(mean(all50(:,17))) '   Right: ' num2str(mean(all50(:,18)))]);
    % display (['ALL      Mean: ' num2str(mean( [all50(:,17) ; all50(:,18)] )) '(' num2str(std( [all50(:,17) ; all50(:,18)] )) ')']);
    %
    % subplot(2,2,1); hold on; title('50cm');
    % hist([all50(:,17) ; all50(:,18)]);
    %
    % display (['100cm means']);
    % display (['rep 1    Left: ' num2str(mean(radials_100_1(:,17))) '   Right: ' num2str(mean(radials_100_1(:,18)))]);
    % display (['rep 2    Left: ' num2str(mean(radials_100_2(:,17))) '   Right: ' num2str(mean(radials_100_2(:,18)))]);
    % display (['rep 3    Left: ' num2str(mean(radials_100_3(:,17))) '   Right: ' num2str(mean(radials_100_3(:,18)))]);
    % display (['rep 4    Left: ' num2str(mean(radials_100_4(:,17))) '   Right: ' num2str(mean(radials_100_4(:,18)))]);
    % display (['ALL      Left: ' num2str(mean(all100(:,17))) '   Right: ' num2str(mean(all100(:,18)))]);
    % display (['ALL      Mean: ' num2str(mean( [all100(:,17) ; all100(:,18)] )) '(' num2str(std( [all100(:,17) ; all100(:,18)] )) ')' ]);
    %
    % subplot(2,2,2); hold on; title('100cm');
    % hist([all100(:,17) ; all100(:,18)]);
    %
    % display (['200cm means']);
    % display (['ALL      Left: ' num2str(mean(radials_200_3(:,17))) '   Right: ' num2str(mean(radials_200_3(:,18)))]);
    % display (['ALL      Mean: ' num2str(mean( [radials_200_3(:,17) ; radials_200_3(:,18)] )) '(' num2str(std( [radials_200_3(:,17) ; radials_200_3(:,18)] )) ')' ]);
    %
    % subplot(2,2,3); hold on; title('200cm');
    % hist([radials_200_3(:,17) ; radials_200_3(:,18)]);
    %
    % display (['450cm means']);
    % display (['rep 1    Left: ' num2str(mean(radials_450_1(:,17))) '   Right: ' num2str(mean(radials_450_1(:,18)))]);
    % display (['rep 2    Left: ' num2str(mean(radials_450_2(:,17))) '   Right: ' num2str(mean(radials_450_2(:,18)))]);
    % display (['rep 3    Left: ' num2str(mean(radials_450_3(:,17))) '   Right: ' num2str(mean(radials_450_3(:,18)))]);
    % display (['rep 4    Left: ' num2str(mean(radials_450_4(:,17))) '   Right: ' num2str(mean(radials_450_4(:,18)))]);
    % display (['ALL      Left: ' num2str(mean(all450(:,17))) '   Right: ' num2str(mean(all450(:,18)))]);
    % display (['ALL      Mean: ' num2str(mean( [all450(:,17) ; all450(:,18)] )) '(' num2str(std( [all450(:,17) ; all450(:,18)] )) ')' ]);
    %
    % subplot(2,2,4); hold on; title('450cm');
    % hist([all450(:,17) ; all450(:,18)]);
    %
    % display (['all target means']);
    % display (['ALL      Left: ' num2str(mean(alltargets(:,17))) '   Right: ' num2str(mean(alltargets(:,18)))]);
    % display (['ALL      Mean: ' num2str(mean( [alltargets(:,17) ; alltargets(:,18)] )) '(' num2str(std( [alltargets(:,17) ; alltargets(:,18)] )) ')' ]);
    %
    % figure(); hold on;
    % subplot(1,2,1); hold on; title('targets');
    % hist([alltargets(:,17) ; alltargets(:,18)]);
    %
    % display (['walk means']);
    % display (['rep 1    Left: ' num2str(mean(walk_1(:,17))) '   Right: ' num2str(mean(walk_1(:,18)))]);
    % display (['rep 2    Left: ' num2str(mean(walk_2(:,17))) '   Right: ' num2str(mean(walk_2(:,18)))]);
    % display (['rep 3    Left: ' num2str(mean(walk_3(:,17))) '   Right: ' num2str(mean(walk_3(:,18)))]);
    % display (['ALL      Left: ' num2str(mean(allwalks(:,17))) '   Right: ' num2str(mean(allwalks(:,18)))]);
    % display (['ALL      Mean: ' num2str(mean( [allwalks(:,17) ; allwalks(:,18)] )) '(' num2str(std( [allwalks(:,17) ; allwalks(:,18)] )) ')' ]);
    %
    % subplot(1,2,2); hold on; title('walks');
    % hist([allwalks(:,17) ; allwalks(:,18)]);
    
    %figure(); hold on;
    %plot(repmat(1,length([alltargets(:,17) ; alltargets(:,18)]),1), [alltargets(:,17) ; alltargets(:,18)], 'm.');
    %h(1) = plot(1,mean( [alltargets(:,17) ; alltargets(:,18)] ),'ko');
    %plot([1 ; 1], [mean( [alltargets(:,17) ; alltargets(:,18)] ) + std( [alltargets(:,17) ; alltargets(:,18)] ) ; mean( [alltargets(:,17) ; alltargets(:,18)] ) - std( [alltargets(:,17) ; alltargets(:,18)] )],'k')
    
    %plot(repmat(2,length([allwalks(:,17) ; allwalks(:,18)]),1), [allwalks(:,17) ; allwalks(:,18)], 'm.');
    %h(2) = plot(2,mean( [allwalks(:,17) ; allwalks(:,18)] ),'ro');
    %plot([2 ; 2], [mean( [allwalks(:,17) ; allwalks(:,18)] ) + std( [allwalks(:,17) ; allwalks(:,18)] ) ; mean( [allwalks(:,17) ; allwalks(:,18)] ) - std( [allwalks(:,17) ; allwalks(:,18)] )],'r')
    
    %xlim([0 3]);
    %ylim([0 3000]);
    %ylabel('pupil area')
    %legend(h,'targets','walks');
    
    pupil_size_le = DPTarget.LE(:,3);
    pupil_size_re = DPTarget.RE(:,3);
    pupil_size_be = [pupil_size_le ; pupil_size_re];
    
    %left eye, right eye target and fixation pts
    fixation_le = bsxfun(@minus,[ DPTarget.LE(:,1:2) repmat(href_dist,size(DPTarget.LE(:,1:2),1),1) ],L);
    target_le = bsxfun(@minus,[ DPTarget.LE(:,4:5) repmat(href_dist,size(DPTarget.LE(:,1:2),1),1) ],L);
    
    fixation_re = bsxfun(@minus,[ DPTarget.RE(:,1:2) repmat(href_dist,size(DPTarget.LE(:,1:2),1),1) ],R);
    target_re = bsxfun(@minus,[ DPTarget.RE(:,4:5) repmat(href_dist,size(DPTarget.LE(:,1:2),1),1) ],R);
    
    %normalize
    n = sqrt(sum(fixation_le.^2,2)); n = n(:,ones(1,3));
    fixation_le_norm = fixation_le./n;
    n = sqrt(sum(fixation_re.^2,2)); n = n(:,ones(1,3));
    fixation_re_norm = fixation_re./n;
    n = sqrt(sum(target_le.^2,2)); n = n(:,ones(1,3));
    target_le_norm = target_le./n;
    n = sqrt(sum(target_re.^2,2)); n = n(:,ones(1,3));
    target_re_norm = target_re./n;
    
    
    error_le = acosd( min( 1,dot( fixation_le_norm , target_le_norm , 2) ) );
    error_re = acosd( min( 1,dot( fixation_re_norm , target_re_norm , 2) ) );
    error_be = [error_le ; error_re ];
    
    subplot(3,1,f); hold on;
    plot(pupil_size_be, error_be, 'ko');
 
    plot([allwalks(:,17) ; allwalks(:,18)], repmat(10,length([allwalks(:,17) ; allwalks(:,18)]),1),'go')
    plot([alltargets(:,17) ; alltargets(:,18)], repmat(12,length([alltargets(:,17) ; alltargets(:,18)]),1),'bo')
    
    %plot([min([allwalks(:,17) ; allwalks(:,18)]) ; min([allwalks(:,17) ; allwalks(:,18)])], [0 25],'m');
    %plot([quantile([allwalks(:,17) ; allwalks(:,18)],.01) ; quantile([allwalks(:,17) ; allwalks(:,18)],.01)], [0 25],'c');
    %plot([quantile([allwalks(:,17) ; allwalks(:,18)],.25) ; quantile([allwalks(:,17) ; allwalks(:,18)],.25)], [0 25],'r');
    %plot([median([allwalks(:,17) ; allwalks(:,18)]) ; median([allwalks(:,17) ; allwalks(:,18)])], [0 25],'k');
    %plot([quantile([allwalks(:,17) ; allwalks(:,18)],.75) ; quantile([allwalks(:,17) ; allwalks(:,18)],.75)], [0 25],'r');
    %plot([quantile([allwalks(:,17) ; allwalks(:,18)],.99) ; quantile([allwalks(:,17) ; allwalks(:,18)],.99)], [0 25],'c');
    %plot([max([allwalks(:,17) ; allwalks(:,18)]) ; max([allwalks(:,17) ; allwalks(:,18)])], [0 25],'m');
    title('error by pupil size, both eyes');
    xlabel('pupil area'); ylabel('error in deg');
    ylim([0 35]);
    
    corr1 = corr( [ error_be(pupil_size_be >= min([allwalks(:,17) ; allwalks(:,18)]) & pupil_size_be <= max([allwalks(:,17) ; allwalks(:,18)])) pupil_size_be(pupil_size_be >= min([allwalks(:,17) ; allwalks(:,18)]) & pupil_size_be <= max([allwalks(:,17) ; allwalks(:,18)])) ] );
    corr2 = corr( [ error_be(pupil_size_be >= quantile([allwalks(:,17) ; allwalks(:,18)],.01) & pupil_size_be <= quantile([allwalks(:,17) ; allwalks(:,18)],.99)) pupil_size_be(pupil_size_be >= quantile([allwalks(:,17) ; allwalks(:,18)],.01) & pupil_size_be <= quantile([allwalks(:,17) ; allwalks(:,18)],.99)) ] );
    corr3 = corr( [ error_be(pupil_size_be >= quantile([allwalks(:,17) ; allwalks(:,18)],.25) & pupil_size_be <= quantile([allwalks(:,17) ; allwalks(:,18)],.75)) pupil_size_be(pupil_size_be >= quantile([allwalks(:,17) ; allwalks(:,18)],.25) & pupil_size_be <= quantile([allwalks(:,17) ; allwalks(:,18)],.75)) ] );
    
    
    text(100,30,['task pupil sizes, error variance (correlation)']);
    text(100,25,['max-min: ' num2str(var(error_be(pupil_size_be >= min([allwalks(:,17) ; allwalks(:,18)]) & pupil_size_be <= max([allwalks(:,17) ; allwalks(:,18)])) )) ' (' num2str(corr1(1,2)) ')'] );
    text(100,20,['1st-99th percentiles: ' num2str(var(error_be(pupil_size_be >= quantile([allwalks(:,17) ; allwalks(:,18)],.01) & pupil_size_be <= quantile([allwalks(:,17) ; allwalks(:,18)],.99)) )) ' (' num2str(corr2(1,2)) ')'] );
    text(100,15,['25th-75th percentiles: ' num2str(var(error_be(pupil_size_be >= quantile([allwalks(:,17) ; allwalks(:,18)],.25) & pupil_size_be <= quantile([allwalks(:,17) ; allwalks(:,18)],.75)) )) ' (' num2str(corr3(1,2)) ')'] );
    
    text(100,10,['ave error size difference, walks-targets:' num2str( ( median(error_be(pupil_size_be >= min([allwalks(:,17) ; allwalks(:,18)]) & pupil_size_be <= max([allwalks(:,17) ; allwalks(:,18)]))) - median(error_be(pupil_size_be >= min([alltargets(:,17) ; alltargets(:,18)]) & pupil_size_be <= max([alltargets(:,17) ; alltargets(:,18)])))  ) ) ]);
    %text(100,23,['correlation between max and min: ' num2str(corr1(1,2))] );
    %text(100,18,['correlation between 1st and 99th percentiles: ' num2str(corr2(1,2))] );
    %text(100,13,['correlation between 25th and 75th percentiles: ' num2str(corr3(1,2))] );
end
keyboard

