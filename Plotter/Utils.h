vector<TString> Split(TString s,TString del){
  TObjArray* array=s.Tokenize(del);
  vector<TString> out;
  for(const auto& obj:*array){
    out.push_back(((TObjString*)obj)->String());
  }
  array->Delete();
  return out;
}

TString Dirname(TString s){
  return s(0,s.Last('/'));
}
TString Basename(TString s){
  return s(s.Last('/')+1,s.Length());
}
