from datetime import datetime
# SET THIS TO FALSE FOR TASK1Bi
test_timing = True
class Node(object):

    def __init__(self, val, next_node=None):
        self.value = val
        self.next_node = next_node

    def get_hash(self):
        if self.value == "A":
            return 0
        elif self.value == "C":
            return 1
        elif self.value == "G":
            return 2
        elif self.value == "T":
            return 3

    def get_next(self):
        return self.next_node

    def set_next(self, n):
        self.next_node = n

    def get_data(self):
        return self.value

    def set_data(self, value):
        self.value = value

    def has_next(self):
        if self.get_next() is None:
            return False
        return True


class LinkedList(object):

    def __init__(self, root=None):
        self.root = root
        self.size = 0
        self.hash_window = None

    def get_size(self):
        return self.size

    # calculate hash for binding site
    def get_hash_value(self,node):
        hash_sum = 0
        for i in range(self.size):
            hash_sum += 4**i * node.get_hash()
            node = node.next_node

        return hash_sum

    def push(self, el):
        # stack way of pushing elements, new element pushed becomes root
        new_node = Node(el, self.root)
        # new el becomes root
        self.root = new_node
        self.size += 1

    # search for your binding site
    def search(self,node, binding_hash_val):
        store_head_pointer = node
        store_node_before_head = node
        head_pointer = node
        tail_pointer = node
        # hash window from head_pointer to tail_pointer
        for steps in range(bindingSiteLength- 1):
            if tail_pointer is not None:
                tail_pointer = tail_pointer.next_node
        while tail_pointer is not None:
            check_hash = head_pointer
            hash_sum = 0

            for steps in range(bindingSiteLength ):
                hash_sum += 4**steps * check_hash.get_hash()
                if check_hash.next_node is not None:
                    check_hash = check_hash.next_node

            if hash_sum == binding_hash_val:
                if head_pointer == self.root:
                    return store_head_pointer
                else:
                    return store_node_before_head

            else:
                store_node_before_head = head_pointer
                head_pointer = head_pointer.next_node
                store_head_pointer = head_pointer
                if tail_pointer is not None:
                    tail_pointer = tail_pointer.next_node
                else:
                    return None
        return None



if __name__ == "__main__":
    # inputs
    if test_timing is False:
        m = input().split()
        DNA = input()
    else:
        total_loading_time = 0.0
        total_finding_time = 0.0
        total_insert_time = 0.0
        start_var_loading = datetime.now()
        data = []
        with open("timing_test_data.txt", "r") as f_in:
            for line in f_in:
                line = line.replace('\n','')
                data.append(line.split(' '))
        m = data[0]
        DNA = data[1][0]
    dnaLength = int(m[0])
    restrictionEnzymes = int(m[1])
    OrigDna = LinkedList()
    for cells in DNA[int(dnaLength)::-1]:
        OrigDna.push(cells)

    for number_of_restriction_enzymes in range(restrictionEnzymes):
        if test_timing is False:
            n = input().split()
        else:
            n = data[2]
        # binding site
        bindingSite = LinkedList()
        bindingSiteLength = int(n[0])
        enzymeLength = int(n[1])

        bindingSiteEls = n[2][:bindingSiteLength][::-1]
        # binding site
        for els in bindingSiteEls:
            bindingSite.push(els)

        splitIndex = int(n[3])
        strand = n[4][:enzymeLength][::-1]
        binding_site_hash_val = bindingSite.get_hash_value(bindingSite.root)
        end_var_loading = datetime.now()
        if test_timing:
            total_loading_time +=(end_var_loading - start_var_loading).total_seconds()
            data.remove(n)

        # binding head search
        start_search = datetime.now()
        bindingHead = OrigDna.search(OrigDna.root, binding_site_hash_val)


        while bindingHead is not None:
            print(bindingHead.value)
            if test_timing:
                found_match = datetime.now()
                total_finding_time += (found_match - start_search).total_seconds()
                load_enzyme = datetime.now()
            # every iteration we need a new enzyme
            new_Enzyme = LinkedList()
            for cells in strand:
                new_Enzyme.push(cells)
            if test_timing:
                finish_load_enzyme = datetime.now()
                total_loading_time += (finish_load_enzyme - load_enzyme).total_seconds()
            # if our binding head is the root node then it means we have to create a new node at the beginning of the linked list
            if bindingHead == OrigDna.root and bindingHead.value == bindingSite.root.value:
                start_insert = datetime.now()
                Head_of_DNA = OrigDna.root
                enzyme_to_be_inserted = new_Enzyme.root
                for steps in range(new_Enzyme.size - 1):
                    enzyme_to_be_inserted = enzyme_to_be_inserted.next_node

                enzyme_to_be_inserted.next_node = Head_of_DNA
                After_enzyme = Head_of_DNA
                # new head of DNA is the root enzyme
                OrigDna.root = new_Enzyme.root
                # timings
                if test_timing:
                    end_insert = datetime.now()
                    total_insert_time += (end_insert - start_insert).total_seconds()
                    start_search = datetime.now()
                bindingHead = OrigDna.search(After_enzyme.next_node, binding_site_hash_val)
                # if splitindex is bigger than 0 we have to move it from our head
                if splitIndex > 0:
                    #timings
                    if test_timing:
                        start_insert = datetime.now()
                    # keep track of head
                    new_Head = After_enzyme
                    Splicing_Position = enzyme_to_be_inserted
                    for steps in range(splitIndex):
                        Splicing_Position = Splicing_Position.next_node
                    if Splicing_Position is not None:
                        After_enzyme = Splicing_Position.next_node
                    else:
                        After_enzyme = Splicing_Position
                    enzyme_to_be_inserted = new_Enzyme.root
                    Splicing_Position.next_node = enzyme_to_be_inserted

                    for steps in range(new_Enzyme.size-1):
                        enzyme_to_be_inserted = enzyme_to_be_inserted.next_node

                    enzyme_to_be_inserted.next_node = After_enzyme
                    After_enzyme = enzyme_to_be_inserted
                    OrigDna.root = new_Head
                    if test_timing:
                        end_insert = datetime.now()
                        total_insert_time += (end_insert - start_insert).total_seconds()
                        start_search = datetime.now()
                    if splitIndex != 0:
                        bindingHead = OrigDna.search(After_enzyme, binding_site_hash_val)
                    else:
                        bindingHead = OrigDna.search(After_enzyme.next_node, binding_site_hash_val)

                # else we insert it at the specific index
            else:
                if test_timing:
                    start_insert = datetime.now()
                for steps in range(splitIndex):
                    if bindingHead.next_node is not None:
                        bindingHead = bindingHead.next_node

                After_enzyme = bindingHead.next_node
                Insert_at_index_0 = bindingHead.next_node
                enzyme_to_be_inserted = new_Enzyme.root
                bindingHead.next_node = enzyme_to_be_inserted
                for enzyme_steps in range(new_Enzyme.size - 1):
                    enzyme_to_be_inserted = enzyme_to_be_inserted.next_node
                enzyme_to_be_inserted.next_node = After_enzyme
                # we resume search from our inserted enzyme
                After_enzyme = enzyme_to_be_inserted
                if test_timing:
                    end_insert = datetime.now()
                    total_insert_time += (end_insert - start_insert).total_seconds()
                    start_search = datetime.now()
                if splitIndex != 0:
                    bindingHead = OrigDna.search(After_enzyme.next_node, binding_site_hash_val)
                else:
                    After_enzyme = Insert_at_index_0
                    bindingHead = OrigDna.search(After_enzyme.next_node, binding_site_hash_val)

    if test_timing:
        print("Total loading time %s" % total_loading_time)
        print("Total finding time %s" % total_finding_time)
        print("Total insert time %s"% total_insert_time)
    else:

        PrintEls = OrigDna.root
        while PrintEls is not None:
            print(PrintEls.value, end="")
            PrintEls = PrintEls.next_node

