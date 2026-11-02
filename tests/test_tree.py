from genefamily.tree import all_nodes, leaves, parse_newick


def test_parse_newick_counts_and_support():
    root = parse_newick("(A:0.1,B:0.2)90:0.05;")
    tips = leaves(root)
    assert [t.name for t in tips] == ["A", "B"]
    assert root.support == 90.0
    assert root.children[0].length == 0.1
    assert len(all_nodes(root)) == 3


def test_parse_newick_handles_nested_tree():
    root = parse_newick("((A:1,B:1)95:1,(C:1,D:1)80:1)100:0;")
    assert [t.name for t in leaves(root)] == ["A", "B", "C", "D"]
    assert root.support == 100.0
