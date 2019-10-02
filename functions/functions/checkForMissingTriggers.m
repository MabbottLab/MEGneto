function p = checkForMissingTriggers( dataset , p)
%Check for missing triggers

%dataset = p.subj.subj_ds{ss};

% extract list of events from the MEG dataset
eventslist = ft_read_event(dataset);

% get list of possible event types
uniquetypes = unique({eventslist.type});

% assign index to events (instead of just string)
eventslistnum = cellfun(@(str) find(strcmp(uniquetypes, str),1), {eventslist.type});

%
markerTest =  find(strcmpi(uniquetypes, 'UPPT002'),1);

includeTrig = find(strcmpi(uniquetypes, char(p.trialdef.details.include{1})),1);
includeOnceTrig = find(strcmpi(uniquetypes, char(p.trialdef.details.includeOnce{1})),1);
excludeTrig = find(strcmpi(uniquetypes, char(p.trialdef.details.includeOnce{1})),1);
t0Trig = find(strcmpi(uniquetypes, p.trialdef.t0marker),1);

includeReplace = find(strcmpi(uniquetypes, char(p.trialdef.details.includeReplace{1})),1);
trigerLength = length([eventslist(eventslistnum == markerTest).sample]);

includeReplace2 = find(strcmpi(uniquetypes, char(p.trialdefError.details.includeReplace{1})),1);

%No trigger found
if isempty(includeTrig)
   fprintf('The total number of triggers for %s were: 0 \n', char(p.trialdef.details.include{1}));
   fprintf('Replacing %s with %s \n\n',char(p.trialdef.details.include{1}), char(p.trialdef.details.includeReplace{1}));
   p.trialdef.Error.includeError = 1;
%Less than expected amount of triggers found
elseif length([eventslist(eventslistnum == includeTrig).sample]) < (trigerLength/2)
   fprintf('The total number of triggers for %s were:  %d \n', char(p.trialdef.details.include{1}),length([eventslist(eventslistnum == includeTrig).sample]));
   fprintf('Replacing %s with %s \n\n',char(p.trialdef.details.include{1}), char(p.trialdef.details.includeReplace{1}));
   p.trialdef.Error.includeError = 1;
end
if isempty(includeOnceTrig)
   fprintf('Tigger %s was not found \n\n', char(p.trialdef.details.includeOnce{1}));
   p.trialdef.Error.includeOnceError = 1; 
elseif length([eventslist(eventslistnum == includeOnceTrig).sample]) < 2
   fprintf('Tigger %s was not found \n\n', char(p.trialdef.details.includeOnce{1}));
   p.trialdef.Error.includeOnceError = 1;
end
if ~isempty(excludeTrig)
    fprintf('Number of excluded triggers found: %d \n\n', length([eventslist(eventslistnum == excludeTrig).sample]));
else
    fprintf('Number of excluded triggers found:0 \n\n');
end

if isempty(t0Trig)
    fprintf('The total number of triggers for t0marker = %s were: 0 \n', p.trialdef.t0marker);
    fprintf('Replacing %s with %s \n\n',p.trialdef.t0marker, char(p.trialdef.details.includeReplace{1}));
    fprintf('Size of includeReplace %d',size(includeReplace));
    if isempty(includeReplace) 
        fprintf('The total number of triggers for Replace-t0marker = %s were: 0 \n', char(p.trialdef.details.includeReplace{1}));
        p.trialdef.Error2.t0marker = 1;
    else
        %fprint('Size of includeReplace %d', length(t0Trig));
        if ~isempty(includeReplace2)
            fprintf('Using Replace-t0marker2 = %s  \n', char(p.trialdefError.details.includeReplace{1}));
            p.trialdef.Error2.t0marker = 1;
        end
    end
elseif length([eventslist(eventslistnum == t0Trig).sample]) < (trigerLength/2)
    fprintf('The total number of triggers for t0marker = %s were: %d \n', p.trialdef.t0marker,length([eventslist(eventslistnum == t0Trig).sample]));
    fprintf('Replacing %s with %s \n\n',p.trialdef.t0marker, char(p.trialdef.details.includeReplace{1}));
    if isempty(includeReplace)
        fprintf('The total number of triggers for Replace-t0marker = %s were: 0 \n', char(p.trialdef.details.includeReplace{1}));
        if ~isempty(includeReplace2)
            fprintf('Using Replace-t0marker2 = %s  \n', char(p.trialdefError.details.includeReplace{1}));
            p.trialdef.Error2.t0marker = 1;
        end
    else
        p.trialdef.Error.t0marker = 1;
    end
    
end




