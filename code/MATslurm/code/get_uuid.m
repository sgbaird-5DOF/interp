function uuid = get_uuid()
% GET_UUID  get unique ID (8 characters, mixture of numbers and letters) via java.util.UUID.randomUUID
temp =  java.util.UUID.randomUUID;
uuidtmp = char(temp.toString);
uuid = uuidtmp(1:8);
end