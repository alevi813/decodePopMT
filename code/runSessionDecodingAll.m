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

newDecodingTesters('direction', {'nancy', 'flat'}, 'nonstimulus', 'nonstimulus');
newDecodingTesters('direction', {'leo', 'flat'}, 'nonstimulus', 'nonstimulus');

newDecodingTesters('direction', {'nancy', 'late'}, 'nonstimulus', 'nonstimulus');
newDecodingTesters('direction', {'leo', 'late'}, 'nonstimulus', 'nonstimulus');

newDecodingTesters('direction', {'nancy', 'early'}, 'nonstimulus', 'nonstimulus');
newDecodingTesters('direction', {'leo', 'early'}, 'nonstimulus', 'nonstimulus');

% choice

newDecodingTesters('choice', {'nancy', 'flat'}, 'nonstimulus', 'nonstimulus');
newDecodingTesters('choice', {'leo', 'flat'}, 'nonstimulus', 'nonstimulus');

newDecodingTesters('choice', {'nancy', 'late'}, 'nonstimulus', 'nonstimulus');
newDecodingTesters('choice', {'leo', 'late'}, 'nonstimulus', 'nonstimulus');

newDecodingTesters('choice', {'nancy', 'early'}, 'nonstimulus', 'nonstimulus');
newDecodingTesters('choice', {'leo', 'early'}, 'nonstimulus', 'nonstimulus');

%% stim nonstim
%direction

newDecodingTesters('direction', {'nancy', 'flat'}, 'stimulus', 'nonstimulus');
newDecodingTesters('direction', {'leo', 'flat'}, 'stimulus', 'nonstimulus');

newDecodingTesters('direction', {'nancy', 'late'}, 'stimulus', 'nonstimulus');
newDecodingTesters('direction', {'leo', 'late'}, 'stimulus', 'nonstimulus');

newDecodingTesters('direction', {'nancy', 'early'}, 'stimulus', 'nonstimulus');
newDecodingTesters('direction', {'leo', 'early'}, 'stimulus', 'nonstimulus');

% choice

newDecodingTesters('choice', {'nancy', 'flat'}, 'stimulus', 'nonstimulus');
newDecodingTesters('choice', {'leo', 'flat'}, 'stimulus', 'nonstimulus');

newDecodingTesters('choice', {'nancy', 'late'}, 'stimulus', 'nonstimulus');
newDecodingTesters('choice', {'leo', 'late'}, 'stimulus', 'nonstimulus');

newDecodingTesters('choice', {'nancy', 'early'}, 'stimulus', 'nonstimulus');
newDecodingTesters('choice', {'leo', 'early'}, 'stimulus', 'nonstimulus');

%% stim/stim

% direction

newDecodingTesters('direction', {'nancy', 'flat'}, 'stimulus', 'stimulus');
newDecodingTesters('direction', {'leo', 'flat'}, 'stimulus', 'stimulus');

newDecodingTesters('direction', {'nancy', 'late'}, 'stimulus', 'stimulus');
newDecodingTesters('direction', {'leo', 'late'}, 'stimulus', 'stimulus');

newDecodingTesters('direction', {'nancy', 'early'}, 'stimulus', 'stimulus');
newDecodingTesters('direction', {'leo', 'early'}, 'stimulus', 'stimulus');

% choice

newDecodingTesters('choice', {'nancy', 'flat'}, 'stimulus', 'stimulus');
newDecodingTesters('choice', {'leo', 'flat'}, 'stimulus', 'stimulus');

newDecodingTesters('choice', {'nancy', 'late'}, 'stimulus', 'stimulus');
newDecodingTesters('choice', {'leo', 'late'}, 'stimulus', 'stimulus');

newDecodingTesters('choice', {'nancy', 'early'}, 'stimulus', 'stimulus');
newDecodingTesters('choice', {'leo', 'early'}, 'stimulus', 'stimulus');

%% nonstim/stim

% direction

newDecodingTesters('direction', {'nancy', 'flat'}, 'nonstimulus', 'stimulus');
newDecodingTesters('direction', {'leo', 'flat'}, 'nonstimulus', 'stimulus');

newDecodingTesters('direction', {'nancy', 'late'}, 'nonstimulus', 'stimulus');
newDecodingTesters('direction', {'leo', 'late'}, 'nonstimulus', 'stimulus');

newDecodingTesters('direction', {'nancy', 'early'}, 'nonstimulus', 'stimulus');
newDecodingTesters('direction', {'leo', 'early'}, 'nonstimulus', 'stimulus');

% choice

newDecodingTesters('choice', {'nancy', 'flat'}, 'nonstimulus', 'stimulus');
newDecodingTesters('choice', {'leo', 'flat'}, 'nonstimulus', 'stimulus');

newDecodingTesters('choice', {'nancy', 'late'}, 'nonstimulus', 'stimulus');
newDecodingTesters('choice', {'leo', 'late'}, 'nonstimulus', 'stimulus');

newDecodingTesters('choice', {'nancy', 'early'}, 'nonstimulus', 'stimulus');
newDecodingTesters('choice', {'leo', 'early'}, 'nonstimulus', 'stimulus');