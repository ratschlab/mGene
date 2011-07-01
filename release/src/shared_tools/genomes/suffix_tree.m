classdef suffix_tree
% class suffix tree
%
% instantiate with 'st = suffix_tree(alphabet)'
%
% type 'methods suffix_tree' for overview of methods
% type 'help suffix_tree.<method_name>' to get method specific help

  properties(SetAccess = 'private', GetAccess = 'private')
    size=0;
    nodes={};
    leafs={};
    sigma={};
    hashf=[];
  end

  methods

    function obj = suffix_tree(alphabet)
    % obj = suffix_tree(alpha)
    % constructs a suffix tree for alphabet in cell <alpha>
    % up to now only alphabets out of consecutive coded letters are allowed
    % e.g. [0:9], ['a':'z'], ...
      obj.sigma=alphabet;
      % here should be a perfect hashfunction created fitting to the alphabet
      % however this works at the moment only for consecutive coded letters e.g [0:9], ['a':'z']...
      [tf,obj.hashf]=find_hash(double(alphabet));
    end

    function n = getsize(obj)
    % returns amount of contained suffices
      n=obj.size;
    end

    function d = getalphabet(obj)
    % returns alphabet
      d=obj.sigma;
    end

    function tf = isempty(obj)
    % tf = isempty(obj)
    % returns true, if no suffices are contained otherwise false
      tf=(obj.size==0);
    end

    function obj = insert(obj, word, pos)
    % [] = insert(obj, word, pos)
    % inserts suffix <word> into the suffix tree
    % word must be a string containing only characters of alphabet given on
    % construction time to object
      pointer=1;
      if obj.size==0
          % nodes need 1 more field as pointer to a leaf
          obj.nodes{1} = zeros(1,length(obj.sigma)+1);
      end
      for k=1:length(word)
        next=obj.nodes{pointer}(obj.hashf(double(word(k)))+1);
        if next==0
            obj.nodes{end+1} = zeros(1,length(obj.sigma)+1);
            obj.nodes{pointer}(obj.hashf(double(word(k)))+1)=length(obj.nodes);
            next=length(obj.nodes);
        end
        pointer=next;
      end
      leafptr=obj.nodes{pointer}(length(obj.sigma)+1);
      if leafptr==0
          % need a new leaf
          obj.leafs{end+1} = [pos];
          obj.nodes{pointer}(length(obj.sigma)+1)=length(obj.leafs);
      else
          % extend existing leaf
          obj.leafs{leafptr}= [obj.leafs{leafptr},pos];
      end
      obj.size=obj.size+1;
    end

    function [tf,pos] = contains(obj, word)
    % [tf,pos] = contains(obj, word)
    % returns tf=true and pos if suffix <word> is contained in the tree <obj> otherwise false
    % pos is a vector of all positions where <word> occures
      pos=[];
      if obj.size==0
        tf=false;
        return;
      end
      pointer=1;
      for k=1:length(word)
        next=obj.nodes{pointer}(obj.hashf(double(word(k)))+1);
        if next==0
          tf= false;
          return;
        else
          pointer=next;
        end
      end
      leafptr=obj.nodes{pointer}(length(obj.sigma)+1);
      if leafptr==0
        tf=false;
        return;
      else
        tf=true;
        pos=obj.leafs{leafptr};
        return;
      end;
    end
    
    function num = contains_prfx(obj, prfx)
    % num = contains_prfx(obj, prfx)
    % returns amount of all contained words with prefix <prfx>
      num=0;
      if obj.size==0
        return;
      end
      pointer=1;
      for k=1:length(prfx)
        next=obj.nodes{pointer}(obj.hashf(double(prfx(k)))+1);
        if next==0
          return;
        else
          pointer=next;
        end
      end
      leafptr=obj.nodes{pointer}(length(obj.sigma)+1);
      if leafptr~=0
        num=length(obj.leafs(leafptr));
      end
      stack=obj.nodes{pointer}(1:end-1);
      while ~isempty(stack)
        if stack(1)==0
          stack=stack(2:end);
        else
          pointer=stack(1);
          stack=[obj.nodes{pointer}(1:end-1),stack(2:end)];
          if obj.nodes{pointer}(end)~=0;
            num=num+length(obj.leafs(obj.nodes{pointer}(end)));
          end
        end
      end
    end

  end

end