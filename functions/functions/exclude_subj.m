function [p]=exclude_subj(p,excluding_subj)

    if isempty(excluding_subj)
        for ss=1:length(p.subj.include)
            p.subj.include(ss)=true;
        end
    else
        for ss=1:length(excluding_subj)
            tmp_index = strmatch(excluding_subj{ss},p.subj.ID);
            if ~isempty(tmp_index)
                p.subj.include(tmp_index)= false;
            else
                fprintf(['\n',excluding_subj{ss},' was not found in the p-structure.\n']);
            end
        end
    end
    % save p-structure
    disp('Saving...');
    save(p.paths.p_strct,'p','-mat','-v7');
    disp('Done.');

end
