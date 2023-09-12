function C = Accuracy_map(C,gt)
 C = bestMap(gt,C);
 ACC = length(find(gt == C))/length(gt);
end