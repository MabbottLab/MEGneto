function sendEmail(step, contacts)

msg = strcat('echo "Done running', " ", step, '!" | mail -s "MEGneto Update"'," ", string(contacts));
for m = msg
    unix(m);
end