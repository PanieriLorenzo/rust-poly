use super::{BaseStore, DynStore, UniStore};

impl<T> BaseStore for Vec<T> {
    type Elem = T;
}

impl<T> UniStore for Vec<T> {}

impl<T> DynStore for Vec<T> {}
