% % %% direction
% % 
% sessionDecoding('direction', {'nancy', 'flat'});
% sessionDecoding('direction', {'leo', 'flat'});
% % 
% % sessionDecoding('direction', {'nancy', 'late'});
% % sessionDecoding('direction', {'leo', 'late'});
% % 
% % sessionDecoding('direction', {'nancy', 'early'});
% % sessionDecoding('direction', {'leo', 'early'});
% % 
% % %% choice
% % 
% sessionDecoding('choice', {'nancy', 'flat'});
% sessionDecoding('choice', {'leo', 'flat'});
% % 
% % sessionDecoding('choice', {'nancy', 'late'});
% % sessionDecoding('choice', {'leo', 'late'});
% % 
% % sessionDecoding('choice', {'nancy', 'early'});
% % sessionDecoding('choice', {'leo', 'early'});''
% 
%  
% %% direction
% 
% sessionDecodingBackwards('direction', {'nancy', 'flat'});
% sessionDecodingBackwards('direction', {'leo', 'flat'});
% 
% sessionDecodingBackwards('direction', {'nancy', 'late'});
% sessionDecodingBackwards('direction', {'leo', 'late'});
% 
% sessionDecodingBackwards('direction', {'nancy', 'early'});
% sessionDecodingBackwards('direction', {'leo', 'early'});
% 
% %% choice
% 
% sessionDecodingBackwards('choice', {'nancy', 'flat'});
% sessionDecodingBackwards('choice', {'leo', 'flat'});
% 
% sessionDecodingBackwards('choice', {'nancy', 'late'});
% sessionDecodingBackwards('choice', {'leo', 'late'});
% 
% sessionDecodingBackwards('choice', {'nancy', 'early'});
% sessionDecodingBackwards('choice', {'leo', 'early'});

%% new function with fitWindow/plotWindow args

%nonstim nonstim

%direction

sessionDecoding('direction', {'nancy', 'flat'}, 'nonstimulus', 'nonstimulus');
sessionDecoding('direction', {'leo', 'flat'}, 'nonstimulus', 'nonstimulus');

sessionDecoding('direction', {'nancy', 'late'}, 'nonstimulus', 'nonstimulus');
sessionDecoding('direction', {'leo', 'late'}, 'nonstimulus', 'nonstimulus');

sessionDecoding('direction', {'nancy', 'early'}, 'nonstimulus', 'nonstimulus');
sessionDecoding('direction', {'leo', 'early'}, 'nonstimulus', 'nonstimulus');

% choice

sessionDecoding('choice', {'nancy', 'flat'}, 'nonstimulus', 'nonstimulus');
sessionDecoding('choice', {'leo', 'flat'}, 'nonstimulus', 'nonstimulus');

sessionDecoding('choice', {'nancy', 'late'}, 'nonstimulus', 'nonstimulus');
sessionDecoding('choice', {'leo', 'late'}, 'nonstimulus', 'nonstimulus');

sessionDecoding('choice', {'nancy', 'early'}, 'nonstimulus', 'nonstimulus');
sessionDecoding('choice', {'leo', 'early'}, 'nonstimulus', 'nonstimulus');


%% stim/stim

% direction

sessionDecoding('direction', {'nancy', 'flat'}, 'stimulus', 'stimulus');
sessionDecoding('direction', {'leo', 'flat'}, 'stimulus', 'stimulus');

sessionDecoding('direction', {'nancy', 'late'}, 'stimulus', 'stimulus');
sessionDecoding('direction', {'leo', 'late'}, 'stimulus', 'stimulus');

sessionDecoding('direction', {'nancy', 'early'}, 'stimulus', 'stimulus');
sessionDecoding('direction', {'leo', 'early'}, 'stimulus', 'stimulus');

% choice

sessionDecoding('choice', {'nancy', 'flat'}, 'stimulus', 'stimulus');
sessionDecoding('choice', {'leo', 'flat'}, 'stimulus', 'stimulus');

sessionDecoding('choice', {'nancy', 'late'}, 'stimulus', 'stimulus');
sessionDecoding('choice', {'leo', 'late'}, 'stimulus', 'stimulus');

sessionDecoding('choice', {'nancy', 'early'}, 'stimulus', 'stimulus');
sessionDecoding('choice', {'leo', 'early'}, 'stimulus', 'stimulus');

%% stim nonstim
%direction

sessionDecoding('direction', {'nancy', 'flat'}, 'stimulus', 'nonstimulus');
sessionDecoding('direction', {'leo', 'flat'}, 'stimulus', 'nonstimulus');

sessionDecoding('direction', {'nancy', 'late'}, 'stimulus', 'nonstimulus');
sessionDecoding('direction', {'leo', 'late'}, 'stimulus', 'nonstimulus');

sessionDecoding('direction', {'nancy', 'early'}, 'stimulus', 'nonstimulus');
sessionDecoding('direction', {'leo', 'early'}, 'stimulus', 'nonstimulus');

% choice

sessionDecoding('choice', {'nancy', 'flat'}, 'stimulus', 'nonstimulus');
sessionDecoding('choice', {'leo', 'flat'}, 'stimulus', 'nonstimulus');

sessionDecoding('choice', {'nancy', 'late'}, 'stimulus', 'nonstimulus');
sessionDecoding('choice', {'leo', 'late'}, 'stimulus', 'nonstimulus');

sessionDecoding('choice', {'nancy', 'early'}, 'stimulus', 'nonstimulus');
sessionDecoding('choice', {'leo', 'early'}, 'stimulus', 'nonstimulus');

%% nonstim/stim

% direction

sessionDecoding('direction', {'nancy', 'flat'}, 'nonstimulus', 'stimulus');
sessionDecoding('direction', {'leo', 'flat'}, 'nonstimulus', 'stimulus');

sessionDecoding('direction', {'nancy', 'late'}, 'nonstimulus', 'stimulus');
sessionDecoding('direction', {'leo', 'late'}, 'nonstimulus', 'stimulus');

sessionDecoding('direction', {'nancy', 'early'}, 'nonstimulus', 'stimulus');
sessionDecoding('direction', {'leo', 'early'}, 'nonstimulus', 'stimulus');

% choice

sessionDecoding('choice', {'nancy', 'flat'}, 'nonstimulus', 'stimulus');
sessionDecoding('choice', {'leo', 'flat'}, 'nonstimulus', 'stimulus');

sessionDecoding('choice', {'nancy', 'late'}, 'nonstimulus', 'stimulus');
sessionDecoding('choice', {'leo', 'late'}, 'nonstimulus', 'stimulus');

sessionDecoding('choice', {'nancy', 'early'}, 'nonstimulus', 'stimulus');
sessionDecoding('choice', {'leo', 'early'}, 'nonstimulus', 'stimulus');